#ifndef DUNE_FEM_CGINVERSEOPERATOR_HH
#define DUNE_FEM_CGINVERSEOPERATOR_HH

#include <dune/common/static_assert.hh>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  namespace Fem
  {

    // ConjugateGradientSolver
    // -----------------------

    /** \class ConjugateGradientSolver
     *  \ingroup OEMSolver
     *  \brief   linear solver using the CG algorithm
     *
     *  \param  Operator  type of the operator to invert
     */
    template< class Operator>
    class ConjugateGradientSolver
    {
      typedef ConjugateGradientSolver< Operator> ThisType;

    public:
      //! type of the operators to invert
      typedef Operator OperatorType;
     
      //! field type of the operator's domain vectors
      typedef typename OperatorType::DomainFieldType DomainFieldType;
      //! field type of the operator's range vectors
      typedef typename OperatorType::RangeFieldType RangeFieldType;
      
      //! type of the operator's domain vectors
      typedef typename OperatorType::DomainFunctionType DomainFunctionType;
      //! type of the operator's range vectors
      typedef typename OperatorType::RangeFunctionType RangeFunctionType;

      //! type of the preconditioner, maps from the range of the operator (the dual space) in it's domain
      typedef Fem::Operator< RangeFunctionType, DomainFunctionType > PreconditionerType;


    private:
      dune_static_assert( (Conversion< DomainFunctionType, RangeFunctionType >::sameType),
                          "DomainFunctionType must equal RangeFunctionType." );

    public:
      /** \brief constructor
       *
       *  \param[in]  epsilon        tolerance
       *  \param[in]  maxIterations  maximum number of CG iterations
       *  \param[in]  verbose        verbose output
       */
      ConjugateGradientSolver ( const RangeFieldType &epsilon,
                                unsigned int maxIterations,
                                bool verbose )
      : epsilon_( epsilon ),
        maxIterations_( maxIterations ),
        verbose_( verbose ),
        averageCommTime_( 0.0 ),
        realCount_( 0 )
        
      {}

      /** \brief constructor
       *
       *  \param[in]  epsilon        tolerance
       *  \param[in]  maxIterations  maximum number of CG iterations
       */
      ConjugateGradientSolver ( RangeFieldType epsilon,
                                unsigned int maxIterations )
      : epsilon_( epsilon ),
        maxIterations_( maxIterations ),
        verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
        averageCommTime_( 0.0 ),
        realCount_( 0 )
      {}

    private:
      // prohibit copying
      ConjugateGradientSolver ( const ThisType & );

    public:
      /** \brief solve \f$op( x ) = b\f$
       *
       *  \note The CG algorithm also works for positive semidefinite operators.
       *        In this case, \f$x \cdot v = b \cdot v\f$ for all \f$v\f$ in the
       *        operator's kernel.
       *
       *  \param[in]   op  linear operator to invert (must be symmetic and
       *                   positive definite)
       *  \param[in]   b   right hand side
       *  \param       x   solution (must be initialized to a start value)
       */
      void solve ( const OperatorType &op, const RangeFunctionType &b, DomainFunctionType &x ) const;

      /** \brief solve \f$op( x ) = b\f$
       *
       *  \note The CG algorithm also works for positive semidefinite operators.
       *        In this case, \f$x \cdot v = b \cdot v\f$ for all \f$v\f$ in the
       *        operator's kernel.
       *
       *  \param[in]   op  linear operator to invert (must be symmetic and
       *                   positive definite)
       *  \param[in]   p   (lef) preconditioning operator  
       *  \param[in]   b   right hand side
       *  \param       x   solution (must be initialized to a start value)
       */
      void solve ( const OperatorType &op, const PreconditionerType &p, 
                   const RangeFunctionType &b, DomainFunctionType &x ) const;

      //! number of iterations needed for last solve 
      unsigned int iterations () const 
      {
        return realCount_;
      }
      
      //! return average communication time during last solve 
      double averageCommTime() const 
      {
        return averageCommTime_;
      }

    protected:
      const RangeFieldType epsilon_;
      const unsigned int maxIterations_;
      const bool verbose_;
      mutable double averageCommTime_;
      mutable unsigned int realCount_;
    };



    // CGInverseOperator
    // -----------------
   
    /** \class   CGInverseOperator
     *  \ingroup OEMSolver
     *  \brief   Inverse operator base on CG method
     */
    template< class DiscreteFunction > 
    class CGInverseOperator
    : public Fem::Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Fem::Operator< DiscreteFunction, DiscreteFunction > BaseType;
      
    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef Fem::Operator< DomainFunctionType, RangeFunctionType > OperatorType;
      typedef Fem::Operator< RangeFunctionType, DomainFunctionType > PreconditionerType;
      
    private:
      typedef ConjugateGradientSolver< OperatorType > SolverType;

    public:
      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       *  \param[in]  verbose  verbosity
       */
      CGInverseOperator ( const OperatorType &op,
                          double redEps, double absLimit,
                          unsigned int maxIter, bool verbose )
        : operator_( op ),
          preconditioner_ ( 0 ),
          solver_( absLimit, maxIter, verbose )
      {}

      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       */
      CGInverseOperator ( const OperatorType &op,
                          double redEps, double absLimit,
                          unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : operator_( op ),
          preconditioner_ ( 0 ),
          solver_( absLimit, maxIter )
      {}
      
      /** \brief constructor of CGInverseOperator
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  redEps   reduction epsilon
       *  \param[in]  absLimit absolut limit of residual
       *  \param[in]  maxIter  maximum number of iteration steps
       */
      CGInverseOperator ( const OperatorType &op,
                          const PreconditionerType &precond,
                          double redEps, double absLimit,
                          unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : operator_( op ),
          preconditioner_( &precond ),
          solver_( absLimit, maxIter )
      {}

      /** \brief application operator
       *
       *  The application operator actually solves the linear system
       *  \f$op(w) = u\f$ using the CG method.
       *
       *  \param[in]   u  argument discrete function
       *  \param[out]  w  destination discrete function
       */
      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        if(preconditioner_)
          solver_.solve( operator_, *preconditioner_, u, w );
        else 
          solver_.solve(operator_,u,w);
      }

      //! number of iterations needed for last solve 
      unsigned int iterations () const 
      {
        return solver_.iterations();
      }
      
      //! return average communication time during last solve 
      double averageCommTime() const 
      {
        return solver_.averageCommTime();
      }

    protected:
      const OperatorType &operator_;
      const PreconditionerType *preconditioner_;
      SolverType solver_;
    };



    // Implementation of ConjugateGradientSolver
    // -----------------------------------------

    template< class Operator >
    inline void ConjugateGradientSolver< Operator >
      ::solve ( const OperatorType &op, const RangeFunctionType &b, DomainFunctionType &x ) const
    {
      const bool verbose = (verbose_ && (b.space().grid().comm().rank() == 0));
      
      const RangeFieldType tolerance = (epsilon_ * epsilon_) * b.scalarProductDofs( b ); 

      averageCommTime_ = 0.0;
      
      RangeFunctionType h( b );
      op( x, h );

      RangeFunctionType r( h );
      r -= b;

      RangeFunctionType p( b );
      p -= h;

      RangeFieldType prevResiduum = 0;
      RangeFieldType residuum = r.scalarProductDofs( r );
   
      for( realCount_ = 0; (residuum > tolerance) && (realCount_ < maxIterations_); ++realCount_ )
      {
        if( realCount_ > 0 )
        { 
          p *= (residuum / prevResiduum);
          p -= r;
        }

        op( p, h );

        const RangeFieldType alpha = residuum / p.scalarProductDofs( h );
        x.axpy( alpha, p );
        r.axpy( alpha, h );

        prevResiduum = residuum;
        residuum = r.scalarProductDofs( r );
        
        double exchangeTime = h.space().communicator().exchangeTime();
        if( verbose )
        {
          std::cerr << "CG-Iteration: " << realCount_ << ", Residuum: " << residuum << std::endl;
          // only for parallel apps 
          if( b.space().grid().comm().size() > 1 )
            std::cerr << "Communication needed: " << exchangeTime << " s" << std::endl;
        }
        
        averageCommTime_ += exchangeTime;
      }
    }


    template< class Operator >
    inline void ConjugateGradientSolver< Operator >
    ::solve ( const OperatorType &op, const PreconditionerType &precond, const RangeFunctionType &b, DomainFunctionType &x ) const
    {
      const bool verbose = (verbose_ && (b.space().grid().comm().rank() == 0));
      
      const RangeFieldType tolerance = (epsilon_ * epsilon_) * b.scalarProductDofs( b ); 

      averageCommTime_ = 0.0;
      
      RangeFunctionType h( b );
      //h=Ax
      op( x, h );

      //r=Ax-b
      RangeFunctionType r( h );
      r -= b;
   
      //p=b-A*x <= r_0 Deufelhard
      RangeFunctionType p( b );
      p -= h;
      
      //q=B*p <=q Deuf
      RangeFunctionType q ( b );
      precond(p,q);
   
      RangeFunctionType s (q);
      
      RangeFieldType prevResiduum = 0;
      RangeFieldType residuum = p.scalarProductDofs( q );//<p,Bp>
   
      for( realCount_ = 0; (residuum > tolerance) && (realCount_ < maxIterations_); ++realCount_ )
      {
        if( realCount_ > 0 )
         { 
           const RangeFieldType beta=residuum/prevResiduum;
           q*=beta;
           q+=(s);
         }
        
        op( q, h );
        
        const RangeFieldType alpha = residuum / q.scalarProductDofs( h );//<p,Bp>/<q,Aq>
        x.axpy( alpha, q );
        
        p.axpy( -alpha, h );//r_k
        
        precond(p,s); //B*r_k
        
        prevResiduum = residuum;//<rk-1,B*rk-1>
        
        residuum = p.scalarProductDofs( s );//<rk,B*rk>
        
        double exchangeTime = h.space().communicator().exchangeTime();
        if( verbose )
        {
          std::cerr << "CG-Iteration: " << realCount_ << ", Residuum: " << residuum << std::endl;
          // only for parallel apps 
          if( b.space().grid().comm().size() > 1 )
            std::cerr << "Communication needed: " << exchangeTime << " s" << std::endl;
        }
        
        averageCommTime_ += exchangeTime;
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_CGINVERSEOPERATOR_HH