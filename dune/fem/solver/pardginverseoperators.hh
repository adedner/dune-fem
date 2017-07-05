#ifndef DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH
#define DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/solver/pardg.hh>
#include <dune/fem/io/parameter.hh>

#ifdef USE_PARDG_ODE_SOLVER

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class DomainFunction, class RangeFunction = DomainFunction >
    class ParDGOperator;



    // ParDGOperator for AdaptiveDiscreteFunction
    // ------------------------------------------

    template< class DomainFunctionSpace, class RangeFunctionSpace >
    class ParDGOperator< AdaptiveDiscreteFunction< DomainFunctionSpace >, AdaptiveDiscreteFunction< RangeFunctionSpace > >
    : public PARDG::Function
    {
      typedef ParDGOperator< AdaptiveDiscreteFunction< DomainFunctionSpace >, AdaptiveDiscreteFunction< RangeFunctionSpace > > ThisType;

      typedef AdaptiveDiscreteFunction< DomainFunctionSpace > DomainFunctionType;
      typedef AdaptiveDiscreteFunction< RangeFunctionSpace > RangeFunctionType;

      typedef typename RangeFunctionSpace :: RangeFieldType  RangeFieldType;

    public:
      typedef Operator< DomainFunctionType, RangeFunctionType > OperatorType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;

      ParDGOperator ( const OperatorType &op, const DomainFunctionSpaceType &domainSpace, const RangeFunctionSpaceType &rangeSpace )
      : operator_( op ),
        domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace )
      {}

      void operator() ( const double *u, double *w, int i = 0 )
      {
        DomainFunctionType uFunction( "ParDGOperator u", domainSpace_, (const RangeFieldType *) u );
        RangeFunctionType  wFunction( "ParDGOperator w", rangeSpace_ , (RangeFieldType *) w );
        operator_( uFunction, wFunction );
      }

      int dim_of_argument( int i = 0 ) const
      {
        assert( i == 0 );
        return domainSpace_.size();
      }

      int dim_of_value ( int i = 0 ) const
      {
        assert( i == 0 );
        return rangeSpace_.size();
      }

    private:
      const OperatorType &operator_;
      const DomainFunctionSpaceType &domainSpace_;
      const RangeFunctionSpaceType &rangeSpace_;
    };




    // ParDGGeneralizedMinResInverseOperator
    // -------------------------------------

    template< class DiscreteFunction >
    class ParDGGeneralizedMinResInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

      typedef ParDGOperator< DiscreteFunction, DiscreteFunction > ParDGOperatorType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Operator< DiscreteFunction, DiscreteFunction > PreconditionerType;

      bool getVerbosity( const ParameterReader& parameter ) const
      {
        return parameter.getValue< bool >( "fem.solver.verbose", false );
      }

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op,
                                              double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                              const ParameterReader &parameter = Parameter::container() )
      : ParDGGeneralizedMinResInverseOperator( op, nullptr, redEps, absLimit, maxIterations, verbose, parameter ) {}

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                                              const ParameterReader &parameter = Parameter::container() )
      : ParDGGeneralizedMinResInverseOperator( op, nullptr, redEps, absLimit, std::numeric_limits< unsigned int >::max(),
                                               getVerbosity( parameter ), parameter )
      {}

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op, double redEps, double absLimit,
                                              unsigned int maxIterations,
                                              const ParameterReader &parameter = Parameter::container() )
      : ParDGGeneralizedMinResInverseOperator( op, nullptr, redEps, absLimit, maxIterations,
                                               getVerbosity( parameter ), parameter ) {}


      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                              double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                              const ParameterReader &parameter = Parameter::container() )
      : ParDGGeneralizedMinResInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, verbose, parameter ) {}

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                              double redEps, double absLimit,
                                              const ParameterReader &parameter = Parameter::container() )
      : ParDGGeneralizedMinResInverseOperator( op, &preconditioner, redEps, absLimit, std::numeric_limits< unsigned int >::max(),
                                               getVerbosity( parameter ), parameter )
      {}

      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                              double redEps, double absLimit, unsigned int maxIterations,
                                              const ParameterReader &parameter = Parameter::container() )
      : ParDGGeneralizedMinResInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations,
                                               getVerbosity( parameter ), parameter ) {}


      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        if( ! std::is_same< typename DiscreteFunction::RangeFieldType, double > ::value )
        {
          DUNE_THROW(Dune::NotImplemented,"ParDGGeneralizedMinResInverseOperator only works for double as RangeFieldType" );
        }

        ParDGOperatorType parDGOperator( operator_, w.space(), u.space() );

        const double* U = (const double *) u.leakPointer() ;
        double* W       = (double * ) w.leakPointer();

        if( preconditioner_ )
        {
          ParDGOperatorType parDGPreconditioner( *preconditioner_, w.space(), w.space() );
          solver_.set_preconditioner( parDGPreconditioner );
          solver_.solve( parDGOperator, W, U );
          solver_.unset_preconditioner();
        }
        else
        {
          solver_.solve( parDGOperator, W, U );
        }
      }

      unsigned int iterations () const
      {
        return solver_.number_of_iterations();
      }

    private:
      ParDGGeneralizedMinResInverseOperator ( const OperatorType &op,
                                              const PreconditionerType *preconditioner,
                                              double redEps, double absLimit,
                                              unsigned int maxIterations, bool verbose,
                                              const ParameterReader &parameter = Parameter::container() )
      : solver_( PARDG::Communicator::instance(), parameter.getValue< int >( "fem.solver.gmres.restart", 20 ) ),
        operator_( op ),
        preconditioner_( preconditioner )
      {
        PARDG::set_tolerance( parameter, solver_, redEps, absLimit, "fem.solver.errormeasure" );

        maxIterations = std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations );
        solver_.set_max_number_of_iterations( int( maxIterations ) );

        // only set output when general verbose mode is enabled
        // (basically to avoid output on every rank)
        if( verbose && Parameter :: verbose() )
        {
          solver_.IterativeSolver::set_output( std::cout );
          solver_.DynamicalObject::set_output( std::cout );
        }
      }

      mutable PARDG::GMRES solver_;
      const OperatorType &operator_;
      const PreconditionerType *preconditioner_;
    };

    // ParDGBiCGStabInverseOperator
    // ----------------------------

    template< class DiscreteFunction >
    class ParDGBiCGStabInverseOperator
    : public Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Operator< DiscreteFunction, DiscreteFunction > BaseType;

      typedef ParDGOperator< DiscreteFunction, DiscreteFunction > ParDGOperatorType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef Operator< DiscreteFunction, DiscreteFunction > PreconditionerType;

      bool getVerbosity( const ParameterReader& parameter ) const
      {
        return parameter.getValue< bool >( "fem.solver.verbose", false );
      }

      ParDGBiCGStabInverseOperator ( const OperatorType &op,
                                     double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                     const ParameterReader &parameter = Parameter::container() )
      : ParDGBiCGStabInverseOperator( op, nullptr, redEps, absLimit, maxIterations, verbose, parameter ) {}

      ParDGBiCGStabInverseOperator ( const OperatorType &op,
                                     double redEps, double absLimit,
                                     const ParameterReader &parameter = Parameter::container() )
      : ParDGBiCGStabInverseOperator( op, nullptr, redEps, absLimit, std::numeric_limits< unsigned int >::max(),
                                      getVerbosity( parameter ), parameter ) {}

      ParDGBiCGStabInverseOperator ( const OperatorType &op,
                                     double redEps, double absLimit, unsigned int maxIterations,
                                     const ParameterReader &parameter = Parameter::container() )
      : ParDGBiCGStabInverseOperator( op, nullptr, redEps, absLimit, maxIterations,
                                      getVerbosity( parameter ), parameter ) {}


      ParDGBiCGStabInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                     double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                     const ParameterReader &parameter = Parameter::container() )
      : ParDGBiCGStabInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations, verbose, parameter ) {}

      ParDGBiCGStabInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                     double redEps, double absLimit,
                                     const ParameterReader &parameter = Parameter::container() )
      : ParDGBiCGStabInverseOperator( op, &preconditioner, redEps, absLimit, std::numeric_limits< unsigned int >::max(),
                                      getVerbosity( parameter ), parameter ) {}

      ParDGBiCGStabInverseOperator ( const OperatorType &op, const PreconditionerType &preconditioner,
                                     double redEps, double absLimit, unsigned int maxIterations,
                                     const ParameterReader &parameter = Parameter::container() )
      : ParDGBiCGStabInverseOperator( op, &preconditioner, redEps, absLimit, maxIterations,
                                      getVerbosity( parameter ), parameter ) {}

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        ParDGOperatorType parDGOperator( operator_, w.space(), u.space() );
        if( preconditioner_ )
        {
          ParDGOperatorType parDGPreconditioner( *preconditioner_, w.space(), w.space() );
          solver_.set_preconditioner( parDGPreconditioner );
          solver_.solve( parDGOperator, w.leakPointer(), u.leakPointer() );
          solver_.unset_preconditioner();
        }
        else
          solver_.solve( parDGOperator, w.leakPointer(), u.leakPointer() );
      }

      unsigned int iterations () const
      {
        return solver_.number_of_iterations();
      }

    private:
      ParDGBiCGStabInverseOperator ( const OperatorType &op, const PreconditionerType *preconditioner,
                                     double redEps, double absLimit, unsigned int maxIterations, bool verbose,
                                     const ParameterReader &parameter )
      : solver_( PARDG::Communicator::instance() ),
        operator_( op ),
        preconditioner_( preconditioner )
      {
        PARDG::set_tolerance( parameter, solver_,redEps, absLimit, "fem.solver.errormeasure" );

        maxIterations = std::min( (unsigned int)std::numeric_limits< int >::max(), maxIterations );
        solver_.set_max_number_of_iterations( int( maxIterations ) );

        // only set output when general verbose mode is enabled
        // (basically to avoid output on every rank)
        if( verbose && Parameter :: verbose() )
        {
          solver_.IterativeSolver::set_output( std::cout );
          solver_.DynamicalObject::set_output( std::cout );
        }
      }

      mutable PARDG::BICGSTAB solver_;
      const OperatorType &operator_;
      const PreconditionerType *preconditioner_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifdef USE_PARDG_ODE_SOLVER

#endif // #ifndef DUNE_FEM_SOLVER_PARDGINVERSEOPERATORS_HH
