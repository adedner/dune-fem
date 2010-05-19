#ifndef RUNGEKUTTA_ODE_SOLVER_HH
#define RUNGEKUTTA_ODE_SOLVER_HH

// inlcude all used headers before, that they don not appear in DuneODE 

//- system includes 
#include <iostream>
#include <cmath>
#include <vector>
#include <pthread.h>
#include <cassert>
#include <sys/times.h>

#include <dune/common/mpihelper.hh>

//- Dune includes 
#include <dune/fem/operator/common/spaceoperatorif.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>

//- include runge kutta ode solver 
#include <dune/fem/solver/rungekutta.hh>

// include headers of PARDG 
#include "pardg.hh"

namespace DuneODE {
struct ODEParameters
: public LocalParameter< ODEParameters, ODEParameters >
{ 
  ODEParameters() : 
    min_it( Parameter::getValue< int >( "fem.ode.miniterations" , 14 ) ),
    max_it( Parameter::getValue< int >( "fem.ode.maxiterations" , 16 ) ),
    sigma( Parameter::getValue< double >( "fem.ode.cflincrease" , 1.1 ) )
  {
  }
  virtual PARDG::IterativeLinearSolver *linearSolver(PARDG::Communicator & comm) const
  {
    int cycles = Parameter::getValue< int >( "fem.ode.gmrescycles" , 15 );
    PARDG::IterativeLinearSolver* solver = new PARDG::GMRES(comm,cycles);
    double tol = Parameter::getValue< double >( "fem.ode.solver.tolerance" , 1e-6 );
    static const std::string errorTypeTable[]
      = { "absolute", "relative" };
    int errorType = Parameter::getEnum( "fem.ode.solver.errormeassure", errorTypeTable, 0 );
    solver->set_tolerance(tol,(errorType==1));
    int maxIter = Parameter::getValue< int >( "fem.ode.solver.iterations" , 1000 );
    solver->set_max_number_of_iterations(maxIter);
    return solver;
  }
  virtual double tolerance() const
  {
    return Parameter::getValue< double >( "fem.ode.tolerance" , 1e-8 );
  }
  virtual int iterations() const
  {
    return Parameter::getValue< int >( "fem.ode.iterations" , 1000 );
  }
  virtual int verbose() const
  {
    static const std::string verboseTypeTable[]
      = { "none", "cfl", "full" };
    return Parameter::getEnum( "fem.ode.verbose" , verboseTypeTable, 0 );
  }
  virtual bool cflFactor( const PARDG::ODESolver &ode,
                          const PARDG::IterativeLinearSolver &solver,
                          bool converged,
                          double &factor) const
  {
    const int iter = solver.number_of_iterations();
    factor = 1.;
    bool changed = false;
    if (converged) 
    {
      if (iter < min_it) 
      {
        factor = sigma;
        changed = true;
      }
      else if (iter > max_it) 
      {
        factor = (double)max_it/(sigma*(double)iter);
        changed = true;
      }
    }
    else
    {
      factor = 0.5;
      changed = true;
    }
    return changed;
  }
  const int min_it,max_it;
  const double sigma;
};
#ifdef USE_PARDG_ODE_SOLVER
template <class Operator>
class OperatorWrapper : public PARDG::Function 
{
  // type of discrete function 
  typedef typename Operator::DestinationType DestinationType;
 public:
  //! constructor 
  OperatorWrapper(Operator& op) 
    : op_(op)
  {}

  //! apply operator applies space operator and creates temporary
  //! discrete function using the memory from outside 
  void operator()(const double *u, double *f, int i = 0) 
  {
    // set actual time of iteration step
    op_.setTime( this->time() );
    
    // call operator apply 
    op_( u, f );
  }

  //! return size of argument 
  int dim_of_argument(int i = 0) const 
  { 
    assert( i == 0 );
    return op_.size();
  }
  
  //! return size of destination  
  int dim_of_value(int i = 0) const 
  { 
    assert( i == 0 );
    return op_.size();
  }

  //! return reference to real operator 
  const Operator& op() const { return op_; }

protected:
  // operator to call 
  Operator& op_;
};


/**
   @ingroup ODESolver
   @{
 **/

/* \brief Base class for ParDG ode solvers */
template<class DestinationImp>
class ParDGOdeSolverBase : public OdeSolverInterface<DestinationImp>
{
  ParDGOdeSolverBase(const ParDGOdeSolverBase& );
public:
  //! type of destination function 
  typedef DestinationImp DestinationType; 

  //! type of discretization operator 
  typedef ODESpaceOperatorInterface<DestinationType> OperatorType;

protected:
  //! constructor
  ParDGOdeSolverBase(Dune::TimeProviderBase& tp, 
                     const int order ) :
    timeProvider_( tp ),
    comm_(PARDG::Communicator::instance()),
    order_( order ),
    initialized_( false ),
    odeSolver_( 0 )
  {
  }
 
public:  
  // initialize time step size 
  void initialize ( const DestinationType& U0 ) 
  {
    // initialized dt on first call
    if ( ! initialized_ )     
    {
      // create instance of ode solver 
      if( ! odeSolver_ ) odeSolver_ = createOdeSolver();

      // apply operator once 
      spaceOperator().initialize( U0 );
      // set initial time step 
      timeProvider_.provideTimeStepEstimate( spaceOperator().timeStepEstimate() );
      // now it's initialized 
      initialized_ = true;
    }
  }

protected:  
  //! destructor  
  ~ParDGOdeSolverBase() { delete odeSolver_; odeSolver_ = 0; }

  //! return reference to ode solver 
  PARDG::ODESolver& odeSolver() 
  {
    assert( odeSolver_ );
    return *odeSolver_;
  }
  
public:
  void description(std::ostream& out) const 
  {
    out << name() << ", steps: " << order_;
    // would be nice to add info on cfl number, but TimeProviderBase doesn't have this
        //<< ", cfl: " << this->timeProvider_.cfl() << "\\\\" <<std::endl;
  }

protected:  
  //! return name of ode solver 
  virtual std::string name()  const = 0; 

  //! return ode solver object 
  virtual PARDG::ODESolver* createOdeSolver() = 0;

  //! return reference to discretization operator for time step initialization 
  virtual const OperatorType& spaceOperator() const = 0;
  
protected:
  Dune::TimeProviderBase& timeProvider_;
  PARDG::Communicator & comm_;
  const int order_;
  bool initialized_;

private:
  // use method odeSolver for access 
  PARDG::ODESolver* odeSolver_;
};

///////////////////////////////////////////////////////
//
//  --ExplicitOdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class ExplicitOdeSolver : 
  public ParDGOdeSolverBase< DestinationImp >
{
  typedef ParDGOdeSolverBase< DestinationImp > BaseType; 

protected:  
  using BaseType :: order_;
  using BaseType :: comm_;
  using BaseType :: initialized_;
  using BaseType :: timeProvider_;
  using BaseType :: odeSolver ;

public:
  using BaseType::initialize;

  typedef typename BaseType :: OperatorType    OperatorType;
  typedef typename BaseType :: DestinationType DestinationType; 

  //! constructor 
  ExplicitOdeSolver(OperatorType& op, 
                    Dune :: TimeProviderBase &tp, 
                    const int order, bool verbose = false) :
    BaseType(tp, order ),
    expl_( op ),
    verbose_( verbose )
  {}

  //! destructor 
  virtual ~ExplicitOdeSolver() {}
 
  //! solve system 
  void solve( DestinationType& U0 ) 
  {
    // initialize 
    if( ! initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ExplicitOdeSolver wasn't initialized before first call!");
    }
    
    // get dt 
    const double dt = timeProvider_.deltaT();
    
    // should be larger then zero 
    assert( dt > 0.0 );
    
    // get time 
    const double time = timeProvider_.time();

    // get leakPointer 
    double* u = U0.leakPointer();
    
    // call ode solver 
    const bool convergence = odeSolver().step(time, dt , u);

    // set time step estimate of operator 
    timeProvider_.provideTimeStepEstimate( expl_.op().timeStepEstimate() );
    
    assert( convergence );
    if(! convergence ) 
    {
      timeProvider_.invalidateTimeStep();
      std::cerr << "No Convergence of ExplicitOdeSolver! \n";
    }
  }

protected:
  //! return name of ode solver 
  virtual std::string name() const { return "ExplicitOdeSolver"; }

  //! for initialization 
  const OperatorType& spaceOperator() const { return expl_.op(); }

  //! create explicit ode solver 
  PARDG::ODESolver* createOdeSolver() 
  {
    PARDG::ODESolver* odeSolver = 0;
    switch ( order_ ) 
    {
      case 1  : odeSolver = new PARDG::ExplicitEuler(comm_, expl_); break;
      case 2  : odeSolver = new PARDG::ExplicitTVD2 (comm_, expl_); break;
      case 3  : odeSolver = new PARDG::ExplicitTVD3 (comm_, expl_); break;
      case 4  : odeSolver = new PARDG::ExplicitRK4  (comm_, expl_); break;
      default : odeSolver = new PARDG::ExplicitBulirschStoer(comm_, expl_, 7);
                std::cerr << "Runge-Kutta method of this order not implemented.\n" 
                          << "Using 7-stage Bulirsch-Stoer scheme.\n"
                          << std::endl;
    }

    if( verbose_ )
    {
      odeSolver->DynamicalObject::set_output(cout);
    }
    assert( odeSolver );
    return odeSolver;
  }
 
protected:
  OperatorWrapper<OperatorType> expl_;
  const bool verbose_;
};

///////////////////////////////////////////////////////
//
//  --ImplicitOdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class ImplicitOdeSolver : 
  public ParDGOdeSolverBase< DestinationImp >
{
  typedef ParDGOdeSolverBase< DestinationImp > BaseType; 
protected:
  using BaseType :: order_;
  using BaseType :: comm_;
  using BaseType :: initialized_;
  using BaseType :: timeProvider_;
  using BaseType :: odeSolver ;

public:
  using BaseType::initialize;

  typedef typename BaseType :: OperatorType    OperatorType;
  typedef typename BaseType :: DestinationType DestinationType; 

  ImplicitOdeSolver(OperatorType& op, 
                    Dune::TimeProviderBase& tp,
                    const int order, 
                    bool verbose ) DUNE_DEPRECATED :
    BaseType(tp, order),
    impl_( op ),
    linsolver_( 0 ),
    param_( ODEParameters().clone() ),
    verbose_( param_->verbose() ),
    cfl_(1.0)
  {
  }

  ImplicitOdeSolver(OperatorType& op, 
                    Dune::TimeProviderBase& tp,
                    const int order,
                    const ODEParameters& parameter= ODEParameters() ) :
    BaseType(tp, order),
    impl_( op ),
    linsolver_( 0 ),
    param_( parameter.clone() ),
    verbose_( parameter.verbose() ),
    cfl_(1.0)
  {
  }

protected:  
  //! return name of ode solver 
  virtual std::string name() const { return "ImplicitOdeSolver"; }

  //! return reference to parameter class 
  const ODEParameters& parameter() const { assert( param_ ); return *param_; }

  //! for initialization 
  virtual const OperatorType& spaceOperator() const { return impl_.op(); }

  //! create implicit ode solver     
  PARDG::ODESolver* createOdeSolver()
  {
    linsolver_ = parameter().linearSolver( comm_ );
    assert( linsolver_ );

    PARDG::DIRK* odeSolver = 0;
    switch ( order_ ) 
    {
      case 1: odeSolver = new PARDG::ImplicitEuler(comm_, impl_); break;
      case 2: odeSolver = new PARDG::Gauss2(comm_, impl_); break;
      case 3: odeSolver = new PARDG::DIRK3 (comm_, impl_); break;
      default : std::cerr << "Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }

    assert( odeSolver );

    odeSolver->set_linear_solver( *linsolver_ );
    odeSolver->set_tolerance( parameter().tolerance() );
    odeSolver->set_max_number_of_iterations( parameter().iterations() );
    
    if( verbose_ == 2 ) 
    {
      odeSolver->IterativeSolver::set_output(cout);
      odeSolver->DynamicalObject::set_output(cout);
    }
    return odeSolver;
  }
  
  virtual ~ImplicitOdeSolver() 
  {
    delete linsolver_; linsolver_ = 0;
    delete param_;     param_ = 0;
  }
  
  virtual int numberOfIterations()
  {
    return static_cast<PARDG::DIRK&> (odeSolver()).number_of_iterations();
  }
  
public:  
  //! solve 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ImplicitOdeSolver wasn't initialized before first call!");
    }

    const double dt   = timeProvider_.deltaT();
    assert( dt > 0.0 );
    const double time = timeProvider_.time();

    // get pointer to solution
    double* u = U0.leakPointer();
      
    const bool convergence = odeSolver().step(time , dt , u);

    double factor;
    bool changed = parameter().cflFactor(odeSolver(),
                                         *(linsolver_),
                                         convergence, factor);
    cfl_ *= factor;
    if (convergence)
    {
      timeProvider_.provideTimeStepEstimate( cfl_ * spaceOperator().timeStepEstimate() );

      if( changed && verbose_ >= 1 )
        derr << " New cfl number is: "<< cfl_ << " (number of iterations ("
             << "linear: " << linsolver_->number_of_iterations() 
             << ", ode: " << numberOfIterations() 
             << ")"
             << std::endl;
    } 
    else 
    {
      timeProvider_.provideTimeStepEstimate( cfl_ * dt );
      timeProvider_.invalidateTimeStep();
     
      // output only in verbose mode 
      if( verbose_ >= 1 )
      {
        derr << "No convergence: New cfl number is "<< cfl_ << std :: endl;
      }
    }

    linsolver_->reset_number_of_iterations();
  }

protected:  
  OperatorWrapper<OperatorType> impl_;
  PARDG::IterativeLinearSolver* linsolver_;
  const ODEParameters* param_;
  const int verbose_;
  double cfl_;
}; // end ImplicitOdeSolver


///////////////////////////////////////////////////////
//
//  --SemiImplicitOdeSolver 
//
///////////////////////////////////////////////////////
template<class DestinationImp>
class SemiImplicitOdeSolver : 
  public ImplicitOdeSolver< DestinationImp >
{
  typedef ImplicitOdeSolver< DestinationImp > BaseType; 
protected:
  using BaseType :: order_;
  using BaseType :: comm_;
  using BaseType :: initialized_;
  using BaseType :: timeProvider_;
  using BaseType :: odeSolver ;
  using BaseType :: linsolver_ ;
  using BaseType :: parameter ;
  using BaseType :: impl_ ;
  using BaseType :: verbose_ ;

public:
  using BaseType::initialize;

  typedef typename BaseType :: OperatorType    OperatorType;
  typedef typename BaseType :: DestinationType DestinationType; 

  SemiImplicitOdeSolver(OperatorType& explOp, 
                        OperatorType& implOp, 
                        Dune::TimeProviderBase& tp,
                        const int order, const bool verbose ) DUNE_DEPRECATED :
    BaseType( implOp, tp, order, ODEParameters() ),
    expl_( explOp )
  {
  }

  SemiImplicitOdeSolver(OperatorType& explOp, 
                        OperatorType& implOp, 
                        Dune::TimeProviderBase& tp,
                        const int order,
                        const ODEParameters& parameter=ODEParameters()) :
    BaseType( implOp, tp, order, parameter ),
    expl_( explOp )
  {
  }

protected:  
  //! return name of ode solver 
  virtual std::string name() const { return "SemiImplicitOdeSolver"; }

  //! for initialization 
  const OperatorType& spaceOperator() const { return expl_.op(); }

  //! create implicit ode solver     
  PARDG::ODESolver* createOdeSolver()
  {
    linsolver_ = parameter().linearSolver( comm_ );
    assert( linsolver_ );

    PARDG::SIRK* odeSolver = 0;
    switch ( order_ ) 
    {
      case 1: odeSolver = new PARDG::SemiImplicitEuler(comm_, impl_, expl_); break;
      case 2: odeSolver = new PARDG::IMEX_SSP222      (comm_, impl_, expl_); break;
      case 3: odeSolver = new PARDG::SIRK33           (comm_, impl_, expl_); break;
      default : std::cerr << "Runge-Kutta method of this order not implemented" 
                          << std::endl;
                abort();
    }

    assert( odeSolver );

    odeSolver->set_linear_solver(*linsolver_);
    odeSolver->set_tolerance( parameter().tolerance() );
    odeSolver->set_max_number_of_iterations( parameter().iterations() );
    
    if( verbose_ == 2 ) 
    {
      odeSolver->IterativeSolver::set_output(cout);
      odeSolver->DynamicalObject::set_output(cout);
    }
    return odeSolver;
  }

  virtual int numberOfIterations()
  {
    return static_cast<PARDG::SIRK&> (odeSolver()).number_of_iterations();
  }
  
protected:  
  OperatorWrapper<OperatorType> expl_;
}; // end SemiImplicitOdeSolver

#endif // USE_PARDG_ODE_SOLVER
 
/**
 @} 
**/

} // end namespace DuneODE

#undef USE_EXTERNAL_BLAS
#endif
