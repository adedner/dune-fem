#ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define DUNE_FEM_NEWTONINVERSEOPERATOR_HH

#include <cfloat>
#include <iostream>
#include <memory>

#include <dune/fem/solver/parameter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>

namespace Dune
{

  namespace Fem
  {

    // NewtonParameter
    // ---------------

    struct NewtonParameter
      : public Dune::Fem::LocalParameter< NewtonParameter, NewtonParameter >
    {
      protected:

      std::shared_ptr<SolverParameter> baseParam_;
      // key prefix, default is fem.ode.newton. (can be overloaded by user)
      const std::string keyPrefix_;

      ParameterReader parameter_;
      public:

      NewtonParameter( const SolverParameter& baseParameter, const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : baseParam_( baseParameter.clone() ),
          keyPrefix_( keyPrefix ),
          parameter_( parameter )
      {}

      explicit NewtonParameter( const SolverParameter& baseParameter, const ParameterReader &parameter = Parameter::container() )
        : baseParam_( baseParameter.clone() ),
          keyPrefix_( "fem.solver.newton." ),
          parameter_( parameter )
      {}

      NewtonParameter( const ParameterReader &parameter = Parameter::container() )
        : baseParam_( nullptr ),
          keyPrefix_( "fem.solver.newton." ),
          parameter_( parameter )
      {}

      NewtonParameter( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : baseParam_( nullptr ),
          keyPrefix_( keyPrefix ),
          parameter_( parameter )
      {}

      const ParameterReader &parameter () const { return parameter_; }

      virtual double toleranceParameter () const
      {
        return parameter_.getValue< double >( keyPrefix_ + "tolerance", 1e-6 );
      }

      virtual double linAbsTolParameter ( const double &tolerance )  const
      {
        return parameter_.getValue< double >(keyPrefix_ +  "linabstol", tolerance / 8 );
      }

      virtual double linReductionParameter ( const double &tolerance ) const
      {
        return parameter_.getValue< double >( keyPrefix_ + "linreduction", tolerance / 8 );
      }

      virtual bool newtonVerbose () const
      {
        const bool v = baseParam_? baseParam_->verbose() : false;
        return parameter_.getValue< bool >(keyPrefix_ +  "verbose", v );
      }

      virtual bool linearSolverVerbose () const
      {
        const bool v = baseParam_? baseParam_->verbose() : false;
        return parameter_.getValue< bool >( keyPrefix_ + "linear.verbose", v );
      }

      virtual int maxIterationsParameter () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxiterations", std::numeric_limits< int >::max() );
      }

      virtual int maxLinearIterationsParameter () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
      }

    };



    // NewtonInverseOperator
    // ---------------------

    /** \class NewtonInverseOperator
     *  \brief inverse operator based on a newton scheme
     *
     *  \tparam  Op      operator to invert (must be a DifferentiableOperator)
     *  \tparam  LInvOp  linear inverse operator
     *
     *  \note Verbosity of the NewtonInverseOperator is controlled via the
     *        paramter <b>fem.solver.newton.verbose</b>; it defaults to
     *        <b>fem.solver.verbose</b>.
     */
    template< class JacobianOperator, class LInvOp >
    class NewtonInverseOperator
    : public Operator< typename JacobianOperator::RangeFunctionType, typename JacobianOperator::DomainFunctionType >
    {
      typedef NewtonInverseOperator< JacobianOperator, LInvOp > ThisType;
      typedef Operator< typename JacobianOperator::RangeFunctionType, typename JacobianOperator::DomainFunctionType > BaseType;

    public:
      //! type of operator's Jacobian
      typedef JacobianOperator JacobianOperatorType;

      //! type of operator to invert
      typedef DifferentiableOperator< JacobianOperatorType > OperatorType;

      //! type of linear inverse operator
      typedef LInvOp LinearInverseOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef NewtonParameter ParametersType;

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */
      NewtonInverseOperator ( const OperatorType &op, const NewtonParameter &parameter )
      : op_( op ),
        parameters_( parameter.clone() ),
        tolerance_( parameter.toleranceParameter() ),
        linAbsTol_( parameter.linAbsTolParameter( tolerance_ ) ),
        linReduction_( parameter.linReductionParameter( tolerance_ ) ),
        verbose_( parameter.newtonVerbose() && MPIManager::rank () == 0 ),
        linVerbose_( parameter.linearSolverVerbose() ),
        maxIterations_( parameter.maxIterationsParameter() ),
        maxLinearIterations_( parameter.maxLinearIterationsParameter() )
      {}

      explicit NewtonInverseOperator ( const OperatorType &op, const ParameterReader &parameter = Parameter::container() )
      : op_( op ),
        parameters_( new ParametersType( parameter ) ),
        tolerance_( parameters_->toleranceParameter() ),
        linAbsTol_( parameters_->linAbsTolParameter( tolerance_ ) ),
        linReduction_( parameters_->linReductionParameter( tolerance_ ) ),
        verbose_( parameters_->newtonVerbose() && MPIManager::rank () == 0 ),
        linVerbose_( parameters_->linearSolverVerbose() ),
        maxIterations_( parameters_->maxIterationsParameter() ),
        maxLinearIterations_( parameters_->maxLinearIterationsParameter() )
      {}

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  epsilon  tolerance for norm of residual
       */
      NewtonInverseOperator ( const OperatorType &op, const DomainFieldType &epsilon,
                              const NewtonParameter &parameter )
      : op_( op ),
        parameters_( parameter.clone() ),
        tolerance_( epsilon ),
        linAbsTol_( parameter.linAbsTolParameter( tolerance_ ) ),
        linReduction_( parameter.linReductionParameter( tolerance_ ) ),
        verbose_( parameter.newtonVerbose() ),
        linVerbose_( parameter.linearSolverVerbose() ),
        maxIterations_( parameter.maxIterationsParameter() ),
        maxLinearIterations_( parameter.maxLinearIterationsParameter() )
      {}

      NewtonInverseOperator ( const OperatorType &op, const DomainFieldType &epsilon,
                              const ParameterReader &parameter = Parameter::container() )
      : op_( op ),
        parameters_( new ParametersType( parameter ) ),
        tolerance_( epsilon ),
        linAbsTol_( parameters_->linAbsTolParameter( tolerance_ ) ),
        linReduction_( parameters_->linReductionParameter( tolerance_ ) ),
        verbose_( parameters_->newtonVerbose() && MPIManager::rank () == 0 ),
        linVerbose_( parameters_->linearSolverVerbose() ),
        maxIterations_( parameters_->maxIterationsParameter() ),
        maxLinearIterations_( parameters_->maxLinearIterationsParameter() )
      {}

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const;

      int iterations () const { return iterations_; }
      int linearIterations () const { return linearIterations_; }

      bool converged () const
      {
        // check for finite |residual| - this also works for -ffinite-math-only (gcc)
        const bool finite = (delta_ < std::numeric_limits< DomainFieldType >::max());
        return finite && (iterations_ < maxIterations_) && (linearIterations_ < maxLinearIterations_);
      }

    protected:
      // hold pointer to jacobian operator, if memory reallocation is needed, the operator should know how to handle this.
      template< class ... Args>
      JacobianOperatorType& jacobian ( Args && ... args ) const
      {
        if( !jOp_ )
          jOp_.reset( new JacobianOperatorType( std::forward< Args >( args ) ... ) );
        return *jOp_;
      }

    private:
      const OperatorType &op_;
      const std::unique_ptr< NewtonParameter > parameters_;

      const double tolerance_, linAbsTol_, linReduction_;
      const bool verbose_;
      const bool linVerbose_;
      const int maxIterations_;
      const int maxLinearIterations_;

      mutable DomainFieldType delta_;
      mutable int iterations_;
      mutable int linearIterations_;
      mutable std::unique_ptr< JacobianOperatorType > jOp_;
    };



    // Implementation of NewtonInverseOperator
    // ---------------------------------------

    template< class JacobianOperator, class LInvOp >
    inline void NewtonInverseOperator< JacobianOperator, LInvOp >
      ::operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      DomainFunctionType residual( u );
      RangeFunctionType dw( w );
      JacobianOperatorType& jOp = jacobian( "jacobianOperator", dw.space(), u.space() );

      // compute initial residual
      op_( w, residual );
      residual -= u;
      delta_ = std::sqrt( residual.scalarProductDofs( residual ) );

      for( iterations_ = 0, linearIterations_ = 0; converged() && (delta_ > tolerance_); ++iterations_ )
      {
        if( verbose_ )
          std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;

        // evaluate operator's jacobian
        op_.jacobian( w, jOp );

        // David: With this factor, the tolerance of CGInverseOp is the absolute
        //        rather than the relative error
        //        (see also dune-fem/dune/fem/solver/inverseoperators.hh)
        const int remLinearIts = maxLinearIterations_ - linearIterations_;
        const LinearInverseOperatorType jInv( jOp, linReduction_, linAbsTol_, remLinearIts, linVerbose_, parameters_->parameter() );

        dw.clear();
        jInv( residual, dw );
        linearIterations_ += jInv.iterations();
        w -= dw;

        op_( w, residual );
        residual -= u;
        delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
      }
      if( verbose_ )
        std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
