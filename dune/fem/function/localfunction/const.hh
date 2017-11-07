#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH

#include <algorithm>
#include <type_traits>
#include <utility>

#include <dune/common/dynvector.hh>

#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declerations
    // -----------------------------

    template< class >
    struct DiscreteFunctionTraits;

    class HasLocalFunction;
    class IsDiscreteFunction;
    struct BindableFunction;

    // BasicConstLocalFunction
    // -----------------------

    template < class BasisFunctionSet, class LocalDofVector >
    class BasicConstLocalFunction
    : public LocalFunction< BasisFunctionSet, LocalDofVector >
    {
      typedef BasicConstLocalFunction< BasisFunctionSet, LocalDofVector >  ThisType;
      typedef LocalFunction< BasisFunctionSet, LocalDofVector > BaseType;

    public:
      //! type of Dof
      typedef typename BaseType::DofType DofType;

      //! type of Entity
      typedef typename BaseType :: EntityType EntityType;

      //! type of BasisFunctionSet
      typedef typename BaseType :: BasisFunctionSetType BasisFunctionSetType;

      //! type of LocalDofVector
      typedef typename BaseType :: LocalDofVectorType LocalDofVectorType;

      //! type of SizeType
      typedef typename BaseType::SizeType SizeType;

      //! default ctor
      BasicConstLocalFunction () {}

      explicit BasicConstLocalFunction ( const BasisFunctionSetType & basisFunctionSet ) : BaseType( basisFunctionSet ) {}

      explicit BasicConstLocalFunction ( const LocalDofVectorType &localDofVector ) : BaseType( localDofVector ) {}

      BasicConstLocalFunction ( const BasisFunctionSetType &basisFunctionSet, const LocalDofVectorType &localDofVector )
      : BaseType( basisFunctionSet, localDofVector )
      {}

      explicit BasicConstLocalFunction ( LocalDofVectorType &&localDofVector ) : BaseType( localDofVector ) {}

      BasicConstLocalFunction ( const BasisFunctionSetType &basisFunctionSet, LocalDofVectorType &&localDofVector )
      : BaseType( basisFunctionSet, localDofVector )
      {}

      BasicConstLocalFunction ( const BaseType &other ) : BaseType( other ) {}

      BasicConstLocalFunction ( const ThisType &other ) : BaseType( static_cast<const BaseType &>( other ) ) {}
      BasicConstLocalFunction ( ThisType && other ) : BaseType( static_cast<BaseType&&>(other) ) {}

      const DofType &operator[] ( SizeType i ) const { return static_cast< const BaseType & >( *this )[ i ]; }
      const DofType &operator[] ( SizeType i ) { return static_cast< const BaseType & >( *this )[ i ]; }

      using BaseType::localDofVector;

   protected:
      using BaseType::clear;
      using BaseType::assign;
      using BaseType::operator +=;
      using BaseType::operator -=;
      using BaseType::axpy;
    };

    /** \ingroup LocalFunction
        \class ConstLocalDiscreteFunction
        \brief A constant local function carrying values for one entity

        A ConstLocalDiscreteFunction is a LocalFunction which is basically doing the same as the
        LocalFunction of a discrete function. The difference is that the local dofs
        are not kept as references but are copied to a local storage.
        Therefore, this is a const local function and any modification of dofs is not
        allowed.

        \note Local DoF numbers correspond directly to array indices. Hence it
        may be more cache efficient to generate a ConstLocalFunction when only a
        const access to the local function is needed.

        \param DiscreteFunction type of the discrete function, the
                                local function shall belong to
     */
    template< class DiscreteFunction >
    class ConstLocalDiscreteFunction
    : public BasicConstLocalFunction<
      typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DiscreteFunctionSpaceType::BasisFunctionSetType,
      Dune::DynamicVector< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DofType,
        typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::LocalDofVectorAllocatorType
      :: template rebind< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > > ::DofType > ::other > >
    {
      typedef ConstLocalDiscreteFunction< DiscreteFunction > ThisType;
      typedef BasicConstLocalFunction< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DiscreteFunctionSpaceType::BasisFunctionSetType,
              Dune::DynamicVector< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DofType,
              typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > > :: LocalDofVectorAllocatorType
              :: template rebind< typename DiscreteFunctionTraits< std::remove_const_t< DiscreteFunction > >::DofType >::other  > >
          BaseType;

    public:
      typedef std::remove_const_t< DiscreteFunction > DiscreteFunctionType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef DiscreteFunctionType GridFunctionType;

      typedef typename BaseType::DofType DofType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::LocalDofVectorType LocalDofVectorType;
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianRangeType JacobianRangeType;
      typedef typename BaseType::HessianRangeType HessianRangeType;

      /** \brief constructor creating a local function without binding it to an
                 entity

          Creates the local function without initializing the fields depending on
          the current entity.

          \note Before using the local function it must be initilized by
          \code
          localFunction.init( entity );
          \endcode

          \param[in] df discrete function the local function shall belong to
       */
      explicit ConstLocalDiscreteFunction ( const DiscreteFunctionType &df )
      : BaseType( LocalDofVectorType( df.localDofVectorAllocator() ) ),
        discreteFunction_( &df )
      {
      }

#if 0 // DEPRECATE
      //! cast a MutableLocalFunction into this one !!! expensive !!!
      ConstLocalDiscreteFunction ( const typename DiscreteFunctionType::LocalFunctionType &localFunction )
      : BaseType( localFunction.basisFunctionSet(), LocalDofVectorType( localFunction.size(), localFunction.discreteFunction().localDofVectorAllocator() ) ),
        discreteFunction_( &localFunction.discreteFunction() )
      {
        std::copy( localFunction.localDofVector().begin(), localFunction.localDofVector().end(), localDofVector().begin() );
      }
#endif

      /** \brief constructor creating a local function and binding it to an
                 entity

          Creates the local function and initilizes the fields depending on the
          current entity. It is not necessary, though allowed, to call init
          before using the discrete function.

          \note The degrees of freedom are not initialized by this function.

          \param[in] df      discrete function the local function shall
                             belong to
          \param[in] entity  entity for initialize the local function to
       */
      ConstLocalDiscreteFunction ( const DiscreteFunctionType &df, const EntityType &entity )
      : BaseType( df.space().basisFunctionSet( entity ), LocalDofVectorType( df.localDofVectorAllocator() )  ),
        discreteFunction_( &df )
      {
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }

      //! copy constructor
      ConstLocalDiscreteFunction ( const ThisType &other )
      : BaseType( static_cast<const BaseType &>( other ) ),
        discreteFunction_( other.discreteFunction_ )
      {}

      //! move constructor
      ConstLocalDiscreteFunction ( ThisType &&other )
      : BaseType( static_cast< BaseType &&>( other ) ),
        discreteFunction_( other.discreteFunction_ )
      {}

      //! destructor
      ~ConstLocalDiscreteFunction ( )
      { unbind(); }

      using BaseType::localDofVector;

      using BaseType::evaluate;
      using BaseType::jacobian;
      using BaseType::hessian;

      /** \brief evaluate the local function
       *
       *  \param[in]   x    evaluation point in local coordinates
       *  \returns          value of the function in the given point
       */
      template< class Point >
      RangeType evaluate ( const Point &p ) const
      {
        RangeType val;
        evaluate( p, val );
        return val;
      }

      /** \brief evaluate Jacobian of the local function
       *
       *  \note Though the Jacobian is evaluated on the reference element, the
       *        return value is the Jacobian with respect to the actual entity.
       *
       *  \param[in]   x    evaluation point in local coordinates
       *  \returns          Jacobian of the function in the evaluation point
       */
      template< class Point >
      JacobianRangeType jacobian ( const Point &p ) const
      {
        JacobianRangeType jac;
        jacobian( p, jac );
        return jac;
      }

      /** \brief evaluate Hessian of the local function
       *
       *  \note Though the Hessian is evaluated on the reference element, the
       *        return value is the Hessian with respect to the actual entity.
       *
       *  \param[in]   x        evaluation point in local coordinates
       *  \returns              Hessian of the function in the evaluation point
       */
      template< class Point >
      HessianRangeType hessian ( const Point &p ) const
      {
        HessianRangeType h;
        hessian( p, h );
        return h;
      }

#if 0 // DEPRECATED
      /** \copydoc Dune::Fem::LocalFunction :: init */
      void init ( const EntityType &entity )
      {
        BaseType::init( discreteFunction().space().basisFunctionSet( entity ) );
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }
#endif

      void bind ( const EntityType &entity )
      {
        BaseType::init( discreteFunction().space().basisFunctionSet( entity ) );
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }
      void unbind () {}

      const DiscreteFunctionType &discreteFunction() const { return *discreteFunction_; }
      const GridFunctionType &gridFunction() const { return discreteFunction(); }

    protected:
      const DiscreteFunctionType* discreteFunction_;
    };



    // ConstLocalFunction
    // ------------------

    namespace Impl
    {

      template< class GF, class = void >
      struct ConstLocalFunction;

      template< class GF >
      struct ConstLocalFunction< GF, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, GF >::value > >
      {
        typedef ConstLocalDiscreteFunction< GF > Type;
      };

      template< class GF >
      struct ConstLocalFunction< GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value && !std::is_base_of< Fem::IsDiscreteFunction, GF >::value  && std::is_class< typename GF::LocalFunctionType >::value > >
      {
        struct Type
          : public GF::LocalFunctionType
        {
          typedef GF GridFunctionType;
          typedef typename GridFunctionType::LocalFunctionType::EntityType EntityType;

          typedef typename GF::LocalFunctionType::RangeType RangeType;
          typedef typename GF::LocalFunctionType::JacobianRangeType JacobianRangeType;
          typedef typename GF::LocalFunctionType::HessianRangeType HessianRangeType;

          explicit Type ( const GridFunctionType &gridFunction )
            : GridFunctionType::LocalFunctionType( gridFunction ),
              gridFunction_( gridFunction )
          {}
          explicit Type ( const GridFunctionType &gridFunction, const EntityType &entity )
            : GridFunctionType::LocalFunctionType( gridFunction ),
              gridFunction_( gridFunction )
          {
            bind(entity);
          }

          using GF::LocalFunctionType::evaluate;
          using GF::LocalFunctionType::jacobian;
          using GF::LocalFunctionType::hessian;
          using GF::LocalFunctionType::bind;

          //! evaluate local function
          template< class Point >
          RangeType evaluate ( const Point &p ) const
          {
            RangeType val;
            evaluate( p, val );
            return val;
          }

          //! jacobian of local function
          template< class Point >
          JacobianRangeType jacobian ( const Point &p ) const
          {
            JacobianRangeType jac;
            jacobian( p, jac );
            return jac;
          }

          //! hessian of local function
          template< class Point >
          HessianRangeType hessian ( const Point &p ) const
          {
            HessianRangeType h;
            hessian( p, h );
            return h;
          }

          void bind ( const EntityType &entity ) { bind( entity ); }
          void unbind () {}

          const GridFunctionType &gridFunction () const { return gridFunction_; }

        private:
          const GridFunctionType &gridFunction_;
        };
      };

      template< class GF >
      struct ConstLocalFunction< GF, std::enable_if_t< std::is_base_of< Fem::BindableFunction, GF >::value && !std::is_base_of< Fem::IsDiscreteFunction, GF >::value > >
      {
        struct Type
        {
          typedef GF GridFunctionType;
          typedef typename GF::EntityType EntityType;
          typedef typename GF::RangeFieldType RangeFieldType;
          typedef typename GF::RangeType RangeType;
          typedef typename GF::JacobianRangeType JacobianRangeType;
          typedef typename GF::HessianRangeType HessianRangeType;

          explicit Type ( const GridFunctionType &gridFunction )
            :  gridFunction_( gridFunction )
          {}
          explicit Type ( const GridFunctionType &gridFunction, const EntityType &entity )
            :  gridFunction_( gridFunction )
          { bind(entity); }

          template <class Point>
          void evaluate(const Point &x, RangeType &ret) const
          {
            gridFunction().evaluate(x,ret);
          }
          template <class Point>
          void jacobian(const Point &x, JacobianRangeType &ret) const
          {
            gridFunction().jacobian(x,ret);
          }
          template <class Point>
          void hessian(const Point &x, HessianRangeType &ret) const
          {
            gridFunction().hessian(x,ret);
          }
          unsigned int order() const { return gridFunction().order(); }

          //! evaluate local function
          template< class Point >
          RangeType evaluate ( const Point &p ) const
          {
            RangeType val;
            evaluate( p, val );
            return val;
          }

          //! jacobian of local function
          template< class Point >
          JacobianRangeType jacobian ( const Point &p ) const
          {
            JacobianRangeType jac;
            jacobian( p, jac );
            return jac;
          }

          //! hessian of local function
          template< class Point >
          HessianRangeType hessian ( const Point &p ) const
          {
            HessianRangeType h;
            hessian( p, h );
            return h;
          }

          template< class Quadrature, class ... Vectors >
          void evaluateQuadrature ( const Quadrature &quad, Vectors & ... values ) const
          {
            static_assert( sizeof...( Vectors ) > 0, "evaluateQuadrature needs to be called with at least one vector." );
            std::ignore = std::make_tuple( ( evaluateQuadrature( quad, values ), 1 ) ... );
          }

          template< class Quadrature, class Vector >
          auto evaluateQuadrature ( const Quadrature &quad, Vector &v ) const
          -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, RangeType >::value >
          {
            for( const auto qp : quad )
              v[ qp.index() ] = evaluate( qp );
          }

          template< class Quadrature, class Vector >
          auto evaluateQuadrature ( const Quadrature &quad, Vector &v ) const
          -> std::enable_if_t< std::is_same< std::decay_t< decltype(v[ 0 ]) >, JacobianRangeType >::value >
          {
            for( const auto qp : quad )
              v[ qp.index() ] = jacobian( qp );
          }

          void bind ( const EntityType &entity ) { gridFunction().bind( entity ); }
          void unbind () { gridFunction().unbind(); }

          const GridFunctionType &gridFunction () const { return gridFunction_; }

        private:
          GridFunctionType &gridFunction () { return gridFunction_; }
          GridFunctionType gridFunction_;
        };
      };
    } // namespace Impl


    template< class GridFunction >
    using ConstLocalFunction = typename Impl::ConstLocalFunction< GridFunction >::Type;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH
