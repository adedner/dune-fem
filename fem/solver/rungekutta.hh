#ifndef RUNGEKUTTA_SOLVER_HH
#define RUNGEKUTTA_SOLVER_HH

// inlcude all used headers before, that they don not appear in DuneODE 

//- system includes 
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

//- Dune includes 
#include <dune/fem/misc/timeprovider.hh>
#include <dune/fem/operator/common/spaceoperatorif.hh>
namespace DuneODE 
{

using namespace Dune;
using namespace std;

/** @addtogroup ODESolver

    \remarks 
    The interface for ODE solvers is defined by the class
    OdeSolverInterface. The interface for discretization operators
    working with the OdeSolvers is described by the class SpaceOperatorInterface.
 @{
 **/

/** \brief Interface class for ODE Solver. */ 
template <class DestinationImp>
class OdeSolverInterface 
{
protected:
  //! cosntructor 
  OdeSolverInterface () {}    
public:
  //! type of destination 
  typedef DestinationImp DestinationType;

  //! destructor 
  virtual ~OdeSolverInterface () {}
  
  /** \brief initialze solver 
      \param[in] arg argument to apply internal operator once for intiail time step estimate 
  */
  virtual void initialize(const DestinationType& arg) = 0;
  
  /** \brief solve \f$\partial_t u = L(u)\f$ where \f$L\f$ is the internal operator.   
      \param[in] u unknown to solve for 
  */
  virtual void solve(DestinationType& u) = 0;
};

/** \brief Base class for explicit RungeKutta ODE solver. */
template<class Operator>
class ExplRungeKuttaBase 
{
public:
  typedef typename Operator::DestinationType DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
  // typedef typename SpaceType :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 
private:
  std::vector< std::vector<double> > a;
  std::vector<double> b;
  std::vector<double> c;
  std::vector<DestinationType*> Upd;
protected:  
  const int ord_;

public:
  /** \brief constructor 
    \param[in] op Operator \f$L\f$ 
    \param[in] tp TimeProvider 
    \param[in] pord polynomial order 
    \param[in] verbose verbosity 
  */
  ExplRungeKuttaBase(Operator& op, TimeProvider& tp, 
                     int pord, bool verbose = true ) :
    a(0),b(0),c(0), Upd(0),
    ord_(pord),
    op_(op),
    tp_(tp),
    initialized_(false)
  {
    assert(ord_>0);
    a.resize(ord_);
    for (int i=0; i<ord_; ++i)
    {
      a[i].resize(ord_);
    }
    b.resize(ord_); 
    c.resize(ord_); 
    
    switch (ord_) {
    case 4 :
      a[0][0]=0.;     a[0][1]=0.;     a[0][2]=0.;    a[0][3]=0.;
      a[1][0]=1.0;    a[1][1]=0.;     a[1][2]=0.;    a[1][3]=0.;
      a[2][0]=0.25;   a[2][1]=0.25;   a[2][2]=0.;    a[2][3]=0.;
      a[3][0]=1./6.;  a[3][1]=1./6.;  a[3][2]=2./3.; a[3][3]=0.;
      b[0]=1./6.;     b[1]=1./6.;     b[2]=2./3.;    b[3]=0.;
      c[0]=0.;        c[1]=1.0;       c[2]=0.5;      c[3]=1.0;
      break;
    case 3 :
      a[0][0]=0.;     a[0][1]=0.;     a[0][2]=0.;
      a[1][0]=1.0;    a[1][1]=0.;     a[1][2]=0.;
      a[2][0]=0.25;   a[2][1]=0.25;   a[2][2]=0.;
      b[0]=1./6.;     b[1]=1./6.;     b[2]=2./3.;
      c[0]=0.;        c[1]=1;         c[2]=0.5;
      break;
    case 2 :
      a[0][0]=0.;     a[0][1]=0.;
      a[1][0]=1.0;    a[1][1]=0.;
      b[0]=0.5;       b[1]=0.5;
      c[0]=0;         c[1]=1;
      break;
    case 1:
      a[0][0]=0.;
      b[0]=1.;
      c[0]=0.;
      break;
    default : std::cerr << "Runge-Kutta method of this order not implemented" 
                        << std::endl;
              abort();
    }

    // create update memory 
    for (int i=0; i<ord_; ++i)
    {
      Upd.push_back(new DestinationType("URK",op_.space()) );
    }
    Upd.push_back(new DestinationType("Ustep",op_.space()) );
  }

  //! destructor 
  ~ExplRungeKuttaBase()
  {
    for(size_t i=0; i<Upd.size(); ++i) 
      delete Upd[i];
  }

  //! apply operator once to get dt estimate 
  void initialize(const DestinationType& U0)
  {
    if( ! initialized_ ) 
    {
      // Compute Steps
      op_(U0, *(Upd[0]));
      initialized_ = true;
    }
  }
  
  //! solve the system 
  void solve(DestinationType& U0) 
  {
    // time might change 
    tp_.unlock();
    
    // get cfl * timeStepEstimate 
    const double dt = tp_.deltaT();
    // get time 
    const double t = tp_.time();

    // Compute Steps
    op_(U0, *(Upd[0]));
    
    for (int i=1; i<ord_; ++i) 
    {
      (Upd[ord_])->assign(U0);
      for (int j=0; j<i ; ++j) 
      {
        (Upd[ord_])->addScaled(*(Upd[j]),(a[i][j]*dt));
      }

      // set new time 
      tp_.setTime( t + c[i]*dt );

      // apply operator 
      op_(*(Upd[ord_]),*(Upd[i]));
    }

    // Perform Update
    for (int j=0; j<ord_; ++j) 
    {
      U0.addScaled(*(Upd[j]),(b[j]*dt));
    }
    
    // restore global time 
    tp_.lock();
  }

protected:
  // operator to solve for 
  const Operator& op_;
  // time provider 
  TimeProvider& tp_;
  // init flag 
  bool initialized_;
};

/** \brief Exlicit RungeKutta ODE solver that also behaves like a time
    stepper. */
template<class Operator>
class ExplRungeKutta : public TimeProvider , 
                       public ExplRungeKuttaBase<Operator> 
{
  typedef ExplRungeKuttaBase<Operator> BaseType;
public:
  typedef typename Operator :: DestinationType DestinationType;
  typedef typename DestinationType :: DiscreteFunctionSpaceType SpaceType;
  typedef typename SpaceType :: GridType :: Traits :: CollectiveCommunication DuneCommunicatorType; 

public:
  /** \brief constructor 
    \param[in] op Operator \f$L\f$ 
    \param[in] pord polynomial order 
    \param[in] cfl cfl number 
    \param[in] verbose verbosity 
  */
  ExplRungeKutta(Operator& op,int pord,double cfl, bool verbose = true ) :
    TimeProvider(0.0,cfl),
    BaseType(op,*this,pord,verbose),
    tp_(op.space().grid().comm(),*this), 
    savetime_(0.0), savestep_(1)
  {
    op.timeProvider(this);
  }
  
  /** \brief constructor 
    \param[in] op Operator \f$L\f$ 
    \param[in] pord polynomial order 
    \param[in] cfl cfl number 
    \param[in] startTime start time of time stepper  
    \param[in] verbose verbosity 
  */
  ExplRungeKutta(Operator& op,int pord,double cfl, double startTime, bool verbose = true ) :
    TimeProvider(startTime,cfl),
    BaseType(op,*this,pord,verbose),
    tp_(op.space().grid().comm(),*this), 
    savetime_(startTime), savestep_(1)
  {
    op.timeProvider(this);
  }

  void initialize(const DestinationType& U0)
  {
    if(! this->initialized_)
    {
      // initialize 
      BaseType :: initialize(U0);
    
      // global min of dt and reset of dtEstimate 
      this->tp_.syncTimeStep();
    }
  }
    
  double solve(typename Operator::DestinationType& U0) 
  {
    initialize( U0 );
    
    // solve ode 
    BaseType :: solve (U0);
    
    // calls setTime ( t + dt ); 
    this->tp_.augmentTime();
    
    // global min of dt and reset of dtEstimate 
    this->tp_.syncTimeStep();
    
    return this->tp_.time();
  }
  
  void printGrid(int nr, const typename Operator::DestinationType& U) 
  {
    if (time()>=savetime_) {
      printSGrid(time(),savestep_*10+nr,this->op_.space(),U);
      ++savestep_;
      savetime_+=0.001;
    }
  }
  
  void printmyInfo(string filename) const {
    std::ostringstream filestream;
    filestream << filename;
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    ofs << "ExplRungeKutta, steps: " << this->ord_ << "\n\n";
    ofs << "                cfl: " << this->tp_.cfl() << "\\\\\n\n";
    ofs.close();
    this->op_.printmyInfo(filename);
  }

private:
  // TimeProvider with communicator 
  ParallelTimeProvider<DuneCommunicatorType> tp_;
  double savetime_;
  int savestep_;
};

/** \brief Exlicit RungeKutta ODE solver. */
template<class DestinationImp>
class ExplicitRungeKuttaSolver : 
  public OdeSolverInterface<DestinationImp> ,
  public ExplRungeKuttaBase<SpaceOperatorInterface<DestinationImp> >  
{
  typedef DestinationImp DestinationType; 
  typedef SpaceOperatorInterface<DestinationImp> OperatorType;
  typedef ExplRungeKuttaBase<OperatorType> BaseType;
 public:
  /** \brief constructor 
    \param[in] op Operator \f$L\f$ 
    \param[in] tp TimeProvider 
    \param[in] pord polynomial order 
    \param[in] verbose verbosity 
  */
  ExplicitRungeKuttaSolver(OperatorType& op, TimeProvider& tp, int pord, bool verbose = false) :
    BaseType(op,tp,pord,verbose),
    timeProvider_(tp)
  {
    // CFL upper estimate 
    // double cfl = 0.45 / (2.0 * pord+1);

    // maximal allowed cfl number 
    // tp.provideCflEstimate(cfl); 
    // assert( tp.cfl() <= 1.0 );

    if(verbose) 
    {
      std::cout << "ExplicitRungeKuttaSolver: cfl = " << tp.cfl() << "!\n";
    } 
  }

  //! destructor 
  virtual ~ExplicitRungeKuttaSolver() {}
  
  //! apply operator once to get dt estimate 
  void initialize(const DestinationType& U0)
  {
    BaseType :: initialize(U0);
  }
  
  //! solve system 
  void solve(DestinationType& U0) 
  {
    // initialize 
    if( ! this->initialized_ ) 
    {
      DUNE_THROW(InvalidStateException,"ExplicitRungeKuttaSolver wasn't initialized before first call!");
    }

    // solve ode 
    BaseType :: solve(U0);
  }

private:
  TimeProvider& timeProvider_;
};

/** @} **/
} // end namespace DuneODE 
#endif
