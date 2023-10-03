#ifndef DUNE_FEM_SCHEMES_MASSLUMPING_HH
#define DUNE_FEM_SCHEMES_MASSLUMPING_HH

// fem includes
#include <dune/fem/schemes/galerkin.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/common/bindguard.hh>
#include <dune/fem/quadrature/lumpingquadrature.hh>

namespace Dune
{

  namespace Fem
  {

#if HAVE_PETSC
    //forward declaration of PetscLinearOperator
    template <typename DomainFunction, typename RangeFunction >
    class PetscLinearOperator ;
#endif

    // GalerkinOperator
    // ----------------

    template< class Integrands, class MassIntegrands, class DomainFunction, class RangeFunction = DomainFunction >
    struct MassLumpingOperator
      : public virtual Operator< DomainFunction, RangeFunction >
    {
      typedef DomainFunction DomainFunctionType;
      typedef RangeFunction RangeFunctionType;


      static_assert( std::is_same< typename DomainFunctionType::GridPartType, typename RangeFunctionType::GridPartType >::value, "DomainFunction and RangeFunction must be defined on the same grid part." );

      typedef typename RangeFunctionType::GridPartType GridPartType;
      typedef ThreadIterator< GridPartType >  ThreadIteratorType;

      typedef Impl::GalerkinOperator< Integrands > GalerkinOperatorImplType;

      template <class Space>
      struct LumpingQuadratureSelector
      {
        typedef typename Space :: GridPartType GridPartType;
        typedef CachingLumpingQuadrature< GridPartType, 0 > InteriorQuadratureType;
        static constexpr bool isLumpingQuadrature = true ;
        // not needed, since lumping is only on element terms
        typedef CachingQuadrature< GridPartType, 1, Capabilities::DefaultQuadrature< Space > :: template DefaultQuadratureTraits  > SurfaceQuadratureType;
      };

      // define Galerkin operator with different quadratures
      typedef Impl::GalerkinOperator< MassIntegrands, LumpingQuadratureSelector > MassOperatorImplType;

      typedef typename RangeFunctionType :: DiscreteFunctionSpaceType   DiscreteFunctionSpaceType;

      template< class... Args >
      explicit MassLumpingOperator ( const GridPartType &gridPart,
                                     const Integrands& integrands,
                                     const MassIntegrands& massIntegrands,
                                     Args &&... args )
        : iterators_( gridPart ),
          impl_( gridPart, std::move( integrands ), std::forward< Args >( args )... ),
          mass_( gridPart, std::move( massIntegrands ), std::forward< Args >( args )... ),
          gridSizeInterior_( 0 ),
          communicate_( true )
      {
        if( mass().model().hasSkeleton() )
          DUNE_THROW(InvalidStateException,"MassLumpingOperator: Mass model cannot have a skeleton term");
      }

      void setCommunicate( const bool communicate )
      {
        communicate_ = communicate;
        if( ! communicate_ && Dune::Fem::Parameter::verbose() )
        {
          std::cout << "MassLumpingOperator::setCommunicate: communicate was disabled!" << std::endl;
        }
      }

      void setQuadratureOrders(unsigned int interior, unsigned int surface)
      {
        // only set quadrature orders for Galerkin part, lumping quadratures are fixed
        size_t size = impl_.size();
        for( size_t i=0; i<size; ++i )
          impl_[ i ].setQuadratureOrders(interior,surface);
      }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const final override
      {
        evaluate( u, w );
      }

      template< class GridFunction >
      void operator() ( const GridFunction &u, RangeFunctionType &w ) const
      {
        evaluate( u, w );
      }

      const GridPartType &gridPart () const { return impl().gridPart(); }

      typedef Integrands ModelType;
      typedef Integrands DirichletModelType;
      ModelType &model() const { return impl().model(); }

      const GalerkinOperatorImplType& impl() const { return (*impl_); }
      const MassOperatorImplType& mass() const { return (*mass_); }

      std::size_t gridSizeInterior () const { return gridSizeInterior_; }

    protected:
      // update number of interior elements as sum over threads
      std::size_t gatherGridSizeInterior () const
      {
        std::size_t gridSizeInterior = 0;
        const size_t size = MPIManager::numThreads();
        for( size_t i=0; i<size; ++i )
          gridSizeInterior += impl_[ i ].gridSizeInterior();
        return gridSizeInterior;
      }

      template< class GridFunction >
      void evaluate( const GridFunction &u, RangeFunctionType &w ) const
      {
        iterators_.update();
        w.clear();

        std::shared_mutex mutex;

        auto doEval = [this, &u, &w, &mutex] ()
        {
          // version with locking
          // add galerkin part
          this->impl().evaluate( u, w, this->iterators_, mutex );
          // add mass lumped part
          this->mass().evaluate( u, w, this->iterators_, mutex );
        };

        bool singleThreadModeError = false ;

        try {
          // execute in parallel
          MPIManager :: run ( doEval );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = gatherGridSizeInterior();
        }
        catch ( const SingleThreadModeError& e )
        {
          singleThreadModeError = true;
        }

        // if error occurred, redo the whole evaluation
        if( singleThreadModeError )
        {
          // reset w from previous entries
          w.clear();

          // re-run in single thread mode if previous attempt failed
          // add galerkin part
          impl().evaluate( u, w, iterators_ );
          // add mass lumped part
          mass().evaluate( u, w, iterators_ );

          // update number of interior elements
          gridSizeInterior_ = impl().gridSizeInterior();
        }

        // synchronize data
        if( communicate_ )
          w.communicate();
      }

      // GalerkinOperator implementation (see galerkin.hh)
      mutable ThreadIteratorType iterators_;
      ThreadSafeValue< GalerkinOperatorImplType > impl_;
      ThreadSafeValue< MassOperatorImplType     > mass_;

      mutable std::size_t gridSizeInterior_;
      bool communicate_;
    };



    // DifferentiableGalerkinOperator
    // ------------------------------

    template< class Integrands, class MassIntegrands, class JacobianOperator >
    class MassLumpingDifferentiableOperator
      : public MassLumpingOperator< Integrands, MassIntegrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType >,
        public DifferentiableOperator< JacobianOperator >
    {
     typedef MassLumpingOperator< Integrands, MassIntegrands, typename JacobianOperator::DomainFunctionType, typename JacobianOperator::RangeFunctionType > BaseType;

    public:
      typedef JacobianOperator JacobianOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainDiscreteFunctionSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeDiscreteFunctionSpaceType;

      typedef DiagonalAndNeighborStencil< DomainDiscreteFunctionSpaceType, RangeDiscreteFunctionSpaceType >  DiagonalAndNeighborStencilType;
      typedef DiagonalStencil< DomainDiscreteFunctionSpaceType, RangeDiscreteFunctionSpaceType >             DiagonalStencilType;

      typedef typename BaseType::GridPartType GridPartType;

      typedef typename BaseType :: GalerkinOperatorImplType  GalerkinOperatorImplType;
      typedef typename BaseType :: MassOperatorImplType      MassOperatorImplType;

      template< class... Args >
      explicit MassLumpingDifferentiableOperator ( const DomainDiscreteFunctionSpaceType &dSpace,
                                                   const RangeDiscreteFunctionSpaceType &rSpace,
                                                   Args &&... args )
        : BaseType( rSpace.gridPart(), std::forward< Args >( args )... ),
          dSpace_(dSpace),
          rSpace_(rSpace),
          domainSpaceSequence_(dSpace.sequence()),
          rangeSpaceSequence_(rSpace.sequence()),
          stencilDAN_(), stencilD_()
      {
        if( impl().model().hasSkeleton() )
          stencilDAN_.reset( new DiagonalAndNeighborStencilType( dSpace_, rSpace_ ) );
        else
          stencilD_.reset( new DiagonalStencilType( dSpace_, rSpace_ ) );
      }

      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const final override
      {
        // assemble Jacobian, same as GalerkinOperator
        assemble( u, jOp );
      }

      template< class GridFunction >
      void jacobian ( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        assemble( u, jOp );
      }

      const DomainDiscreteFunctionSpaceType& domainSpace() const
      {
        return dSpace_;
      }
      const RangeDiscreteFunctionSpaceType& rangeSpace() const
      {
        return rSpace_;
      }
      // needed for DGHelmholtz operator
      const RangeDiscreteFunctionSpaceType& space() const
      {
        return rangeSpace();
      }

      using BaseType::gridPart;

    protected:
      using BaseType::impl;
      using BaseType::mass;
      using BaseType::gatherGridSizeInterior;
      using BaseType::iterators_;
      using BaseType::gridSizeInterior_;

      void prepare( JacobianOperatorType& jOp ) const
      {
        if ( domainSpaceSequence_ != domainSpace().sequence()
             || rangeSpaceSequence_ != rangeSpace().sequence() )
        {
          domainSpaceSequence_ = domainSpace().sequence();
          rangeSpaceSequence_ = rangeSpace().sequence();
          if( impl().model().hasSkeleton() )
          {
            assert( stencilDAN_ );
            stencilDAN_->update();
          }
          else
          {
            assert( stencilD_ );
            stencilD_->update();
          }
        }

        if( impl().model().hasSkeleton() )
          jOp.reserve( *stencilDAN_ );
        else
          jOp.reserve( *stencilD_ );
        // set all entries to zero
        jOp.clear();
      }

      template < class GridFunction >
      void assemble( const GridFunction &u, JacobianOperatorType &jOp ) const
      {
        prepare( jOp );
        iterators_.update();

        bool singleThreadModeError = false;
        std::shared_mutex mutex;

        auto doAssemble = [this, &u, &jOp, &mutex] ()
        {
          std::tuple< const GalerkinOperatorImplType&, const MassOperatorImplType& > ops( this->impl(), this->mass() );
          this->impl().assemble( u, jOp, this->iterators_, ops, mutex );
          // assemble Jacobian, same as GalerkinOperator
          //this->impl().assemble( u, jOp, this->iterators_, mutex );

          // add mass lumped terms
          //this->mass().assemble( u, jOp, this->iterators_, mutex );
        };

        try {
          // execute in parallel
          MPIManager :: run ( doAssemble );

          // update number of interior elements as sum over threads
          gridSizeInterior_ = gatherGridSizeInterior();
        }
        catch ( const SingleThreadModeError& e )
        {
          singleThreadModeError = true;
        }

        if (singleThreadModeError)
        {
          // redo matrix assembly since it failed
          jOp.clear();
          std::tuple< const GalerkinOperatorImplType&, const MassOperatorImplType& > ops( this->impl(), this->mass() );
          this->impl().assemble( u, jOp, this->iterators_, ops, mutex );
          // add galerkin terms
          //impl().assemble( u, jOp, iterators_ );
          // add mass lumped terms
          //mass().assemble( u, jOp, this->iterators_, mutex );

          // update number of interior elements
          gridSizeInterior_ = impl().gridSizeInterior();
        }

        // note: assembly done without local contributions so need
        // to call flush assembly
        jOp.flushAssembly();
      }

      const DomainDiscreteFunctionSpaceType &dSpace_;
      const RangeDiscreteFunctionSpaceType &rSpace_;
      mutable int domainSpaceSequence_, rangeSpaceSequence_;

      mutable std::unique_ptr< DiagonalAndNeighborStencilType > stencilDAN_;
      mutable std::unique_ptr< DiagonalStencilType >            stencilD_;
    };



    // AutomaticDifferenceGalerkinOperator
    // -----------------------------------

    template< class Integrands, class DomainFunction, class RangeFunction >
    class MassLumpingAutomaticDifferenceGalerkinOperator
      : public MassLumpingOperator< Integrands, DomainFunction, RangeFunction >,
        public AutomaticDifferenceOperator< DomainFunction, RangeFunction >
    {
      typedef MassLumpingOperator< Integrands, DomainFunction, RangeFunction > BaseType;
      typedef AutomaticDifferenceOperator< DomainFunction, RangeFunction > AutomaticDifferenceOperatorType;

    public:
      typedef typename BaseType::GridPartType GridPartType;

      template< class... Args >
      explicit MassLumpingAutomaticDifferenceGalerkinOperator ( const GridPartType &gridPart, Args &&... args )
        : BaseType( gridPart, std::forward< Args >( args )... ), AutomaticDifferenceOperatorType()
      {}
    };



    namespace Impl
    {

      // MassLumpingSchemeImpl
      // ---------------------
      template< class Integrands, class MassIntegrands,
                class LinearOperator, class InverseOperator, bool addDirichletBC,
                template <class,class,class> class DifferentiableGalerkinOperatorImpl = MassLumpingDifferentiableOperator >
      struct MassLumpingSchemeImpl
      {
        typedef InverseOperator InverseOperatorType;
        typedef Integrands      ModelType;
        typedef MassIntegrands  MassModelType;

        using DifferentiableOperatorType = std::conditional_t< addDirichletBC,
           DirichletWrapperOperator< DifferentiableGalerkinOperatorImpl< Integrands, MassIntegrands, LinearOperator >>,
           DifferentiableGalerkinOperatorImpl< Integrands, MassIntegrands, LinearOperator > >;

        using DirichletBlockVector = typename DirichletBlockSelector<
                 DirichletWrapperOperator<
                    DifferentiableGalerkinOperatorImpl< Integrands, MassIntegrands, LinearOperator >>,
                 addDirichletBC>::type;

        typedef typename DifferentiableOperatorType::DomainFunctionType DomainFunctionType;
        typedef typename DifferentiableOperatorType::RangeFunctionType RangeFunctionType;
        typedef typename DifferentiableOperatorType::JacobianOperatorType LinearOperatorType;
        typedef typename DifferentiableOperatorType::JacobianOperatorType JacobianOperatorType;

        typedef RangeFunctionType DiscreteFunctionType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeFunctionSpaceType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType DomainFunctionSpaceType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

        typedef typename DiscreteFunctionSpaceType::FunctionSpaceType FunctionSpaceType;
        typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

        typedef Dune::Fem::NewtonInverseOperator< LinearOperatorType, InverseOperator > NewtonOperatorType;
        typedef InverseOperator LinearInverseOperatorType;
        typedef typename NewtonOperatorType::ErrorMeasureType ErrorMeasureType;

        typedef PreconditionerFunctionWrapper<
              typename LinearOperatorType::RangeFunctionType,
              typename LinearOperatorType::DomainFunctionType >  PreconditionerFunctionWrapperType;

        // std::function to represents the Python function passed as potential preconditioner
        typedef typename PreconditionerFunctionWrapperType::PreconditionerFunctionType  PreconditionerFunctionType ;

        struct SolverInfo
        {
          SolverInfo ( bool converged, int linearIterations, int nonlinearIterations, const std::vector<double>& timing )
            : converged( converged ), linearIterations( linearIterations ),
              nonlinearIterations( nonlinearIterations ), timing( timing )
          {}

          bool converged;
          int linearIterations, nonlinearIterations;
          std::vector<double> timing;
        };

        MassLumpingSchemeImpl ( const DiscreteFunctionSpaceType &dfSpace,
                                const Integrands &integrands,
                                const MassIntegrands& massIntegrands,
                                const ParameterReader& parameter = Parameter::container() )
          : dfSpace_( dfSpace ),
            fullOperator_( dfSpace, dfSpace, std::move(integrands), std::move(massIntegrands) ),
            invOp_(parameter),
            parameter_(parameter)
        {}

        void setQuadratureOrders(unsigned int interior, unsigned int surface) { fullOperator().setQuadratureOrders(interior,surface); }

        const DifferentiableOperatorType &fullOperator() const { return fullOperator_; }
        DifferentiableOperatorType &fullOperator() { return fullOperator_; }

        void constraint ( DiscreteFunctionType &u ) const {}

        template< class GridFunction >
        void operator() ( const GridFunction &u, DiscreteFunctionType &w ) const
        {
          fullOperator()( u, w );
        }

        void setErrorMeasure(ErrorMeasureType &errorMeasure) const
        {
          invOp_.setErrorMeasure(errorMeasure);
        }

        SolverInfo solve ( const DiscreteFunctionType &rhs, DiscreteFunctionType &solution ) const
        {
          DiscreteFunctionType rhs0 = rhs;
          setZeroConstraints( rhs0 );
          setModelConstraints( solution );

          invOp_.bind(fullOperator());
          invOp_( rhs0, solution );
          invOp_.unbind();
          return SolverInfo( invOp_.converged(), invOp_.linearIterations(), invOp_.iterations(), invOp_.timing() );
        }

        SolverInfo solve ( DiscreteFunctionType &solution ) const
        {
          DiscreteFunctionType bnd( solution );
          bnd.clear();
          setModelConstraints( solution );
          invOp_.bind(fullOperator());
          invOp_( bnd, solution );
          invOp_.unbind();
          return SolverInfo( invOp_.converged(), invOp_.linearIterations(), invOp_.iterations(), invOp_.timing() );
        }

        SolverInfo solve ( DiscreteFunctionType &solution, const PreconditionerFunctionType& p ) const
        {
          DiscreteFunctionType bnd( solution );
          bnd.clear();
          setModelConstraints( solution );

          PreconditionerFunctionWrapperType pre( p );
          invOp_.bind(fullOperator(), pre);
          invOp_( bnd, solution );
          invOp_.unbind();
          return SolverInfo( invOp_.converged(), invOp_.linearIterations(), invOp_.iterations(), invOp_.timing() );
        }

        template< class GridFunction >
        void jacobian( const GridFunction &ubar, LinearOperatorType &linearOp) const
        {
          fullOperator().jacobian( ubar, linearOp );
        }

        const DiscreteFunctionSpaceType &space () const { return dfSpace_; }
        const GridPartType &gridPart () const { return space().gridPart(); }
        ModelType &model() const { return fullOperator().model(); }

        void setConstraints( DomainFunctionType &u ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( u );
        }
        void setConstraints( const typename DiscreteFunctionType::RangeType &value, DiscreteFunctionType &u ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( value, u );
        }
        void setConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( u, v );
        }
        void subConstraints( const DiscreteFunctionType &u, DiscreteFunctionType &v ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().subConstraints( u, v );
        }
        const auto& dirichletBlocks() const
        {
          if constexpr (addDirichletBC)
            return fullOperator().dirichletBlocks();
        }

        const ParameterReader& parameter () const
        {
          return parameter_;
        }

      protected:
        void setZeroConstraints( DiscreteFunctionType &u ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( typename DiscreteFunctionType::RangeType(0), u );
        }
        void setModelConstraints( DiscreteFunctionType &u ) const
        {
          if constexpr (addDirichletBC)
            fullOperator().setConstraints( u );
        }
        const DiscreteFunctionSpaceType &dfSpace_;
        DifferentiableOperatorType fullOperator_;
        mutable NewtonOperatorType invOp_;
        const ParameterReader parameter_;
      };

    } // end namespace Impl

    // MassLumpingScheme
    // -----------------

    template< class Integrands, class MassIntegrands, class LinearOperator, class InverseOperator, bool addDirichletBC >
    using MassLumpingScheme = Impl::MassLumpingSchemeImpl< Integrands, MassIntegrands, LinearOperator, InverseOperator, addDirichletBC,
                                                           MassLumpingDifferentiableOperator >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SCHEMES_MOLGALERKIN_HH
