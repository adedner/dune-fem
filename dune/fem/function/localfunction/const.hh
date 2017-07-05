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

      using BaseType::operator[];
      using BaseType::localDofVector;

   protected:
      DofType &operator[] ( SizeType num )
      {
        return static_cast< BaseType &>( *this )[ num ];
      }

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

      //! cast a MutableLocalFunction into this one !!! expensive !!!
      ConstLocalDiscreteFunction ( const typename DiscreteFunctionType::LocalFunctionType &localFunction )
      : BaseType( localFunction.basisFunctionSet(), LocalDofVectorType( localFunction.size(), localFunction.discreteFunction().localDofVectorAllocator() ) ),
        discreteFunction_( &localFunction.discreteFunction() )
      {
        std::copy( localFunction.localDofVector().begin(), localFunction.localDofVector().end(), localDofVector().begin() );
      }

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

      using BaseType::localDofVector;

      /** \copydoc Dune::Fem::LocalFunction :: init */
      void init ( const EntityType &entity )
      {
        BaseType::init( discreteFunction().space().basisFunctionSet( entity ) );
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }

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
      struct ConstLocalFunction< GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value && std::is_base_of< Fem::IsDiscreteFunction, GF >::value > >
      {
        typedef ConstLocalDiscreteFunction< GF > Type;
      };

      template< class GF >
      struct ConstLocalFunction< GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value && !std::is_base_of< Fem::IsDiscreteFunction, GF >::value > >
      {
        struct Type
          : public GF::LocalFunctionType
        {
          typedef GF GridFunctionType;

          explicit Type ( const GridFunctionType &gridFunction )
            : GridFunctionType::LocalFunctionType( gridFunction ),
              gridFunction_( gridFunction )
          {}

          const GridFunctionType &gridFunction () const { return gridFunction_; }

        private:
          const GridFunctionType &gridFunction_;
        };
      };

    } // namespace Impl


    template< class GridFunction >
    using ConstLocalFunction = typename Impl::ConstLocalFunction< GridFunction >::Type;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH
