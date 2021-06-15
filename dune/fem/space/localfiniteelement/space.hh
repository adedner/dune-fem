#ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_SPACE_HH
#define DUNE_FEM_SPACE_LOCALFINITEELEMENT_SPACE_HH

#include <cassert>

#include <memory>
#include <utility>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/fem/common/hybrid.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/basisfunctionset/transformed.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/localinterpolation.hh>
#include <dune/fem/space/mapper/compile.hh>
#include <dune/fem/space/mapper/indexsetdofmapper.hh>
#include <dune/fem/space/shapefunctionset/proxy.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

#include <dune/fem/space/localfiniteelement/shapefunctionset.hh>
#include <dune/fem/space/localfiniteelement/capabilities.hh>
#include <dune/fem/space/localfiniteelement/interpolation.hh>
#include <dune/fem/space/localfiniteelement/localrestrictprolong.hh>

namespace Dune
{

  namespace Fem
  {

    // LocalFiniteElementSpaceTraits
    // -----------------------------

    template< class LFEMap, class FunctionSpace, class Storage >
    struct LocalFiniteElementSpaceTraits
    {
      typedef LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > DiscreteFunctionSpaceType;

      typedef LFEMap LFEMapType;

      typedef typename LFEMapType::GridPartType GridPartType;
      typedef typename LFEMapType::LocalFiniteElementType LocalFiniteElementType;

      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static constexpr int codimension = 0;
      static constexpr bool isScalar = LocalFiniteElementType::Traits::LocalBasisType::Traits::dimRange==1;

      typedef std::conditional_t<isScalar,
              Hybrid::IndexRange< int, FunctionSpace::dimRange >,
              Hybrid::IndexRange< int, 1 >
              > LocalBlockIndices;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      // -1 is default value if pointSetId not available
      static const int pointSetId = detail::SelectPointSetId< LFEMap >::value;

    public:
      // typedef Dune::Fem::IndexSetDofMapper< GridPartType, LagrangeLocalDofMapping< GridPartType > > BlockMapperType;
      typedef Dune::Fem::IndexSetDofMapper< GridPartType > BlockMapperType;

      typedef LocalFunctionsShapeFunctionSet< typename LocalFiniteElementType::Traits::LocalBasisType, pointSetId > LocalFunctionsShapeFunctionSetType;
      typedef SelectCachingShapeFunctionSet< LocalFunctionsShapeFunctionSetType, Storage > StoredShapeFunctionSetType;

      typedef ShapeFunctionSetProxy< StoredShapeFunctionSetType > ShapeFunctionSetProxyType;
      // only extend to vector valued in case that the original space is scalar
      typedef std::conditional_t<isScalar,
              VectorialShapeFunctionSet< ShapeFunctionSetProxyType, typename FunctionSpaceType::RangeType >,
              ShapeFunctionSetProxyType
              > ShapeFunctionSetType;

    private:
      template< class LFEM >
      static TransformedBasisFunctionSet< EntityType, ShapeFunctionSetType, typename LFEM::TransformationType > basisFunctionSet ( const LFEM & );

      static DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > basisFunctionSet ( ... );

    public:
      typedef decltype( basisFunctionSet( std::declval< const LFEMapType & >() ) ) BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // LocalFiniteElementSpace
    // -----------------------

    /** \brief Rannacher-Turek Space
     *  \class LocalFiniteElementSpace
     *  \ingroup DiscreteFunctionSpace
     *
     *  \note The LocalFiniteElementSpace depends on
     *        dune-localfunctions (see http://www.dune-project.org).
     *
     *  \todo please doc me
     **/
    template< class LFEMap, class FunctionSpace, class Storage >
    class LocalFiniteElementSpace
      : public DiscreteFunctionSpaceDefault< LocalFiniteElementSpaceTraits< LFEMap, FunctionSpace, Storage > >
    {
      typedef LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > ThisType;
      typedef DiscreteFunctionSpaceDefault< LocalFiniteElementSpaceTraits< LFEMap, FunctionSpace, Storage > >
        BaseType;

    public:
      typedef typename BaseType::Traits Traits;
      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::Traits::ShapeFunctionSetType ShapeFunctionSetType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef typename BaseType::BlockMapperType BlockMapperType;

      typedef typename Traits::LocalFiniteElementType LocalFiniteElementType;

      typedef typename Traits::LFEMapType LFEMapType;

    private:
      typedef typename LocalFiniteElementType::Traits::LocalBasisType           LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType   LocalInterpolationType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType    LocalCoefficientsType;

      typedef typename LFEMapType::KeyType KeyType;

      typedef typename Traits::LocalFunctionsShapeFunctionSetType LocalFunctionsShapeFunctionSetType;

      struct LFEMapFactory
      {
        static LFEMapType *createObject ( std::pair< GridPartType *, KeyType > key ) { return new LFEMapType( *key.first, key.second ); }
        static void deleteObject ( LFEMapType *object ) { delete object; }
      };

      typedef SingletonList< std::pair< GridPartType *, KeyType >, LFEMapType, LFEMapFactory > LFEMapProviderType;

      typedef typename Traits::StoredShapeFunctionSetType StoredShapeFunctionSetType;
      typedef std::vector< std::unique_ptr< StoredShapeFunctionSetType > > StoredShapeFunctionSetVectorType;

      struct StoredShapeFunctionSetVectorFactory
      {
        static StoredShapeFunctionSetVectorType *createObject ( LFEMapType *lfeMap ) { return new StoredShapeFunctionSetVectorType( lfeMap->size() ); }
        static void deleteObject ( StoredShapeFunctionSetVectorType *object ) { delete object; }
      };

      typedef SingletonList< LFEMapType *, StoredShapeFunctionSetVectorType, StoredShapeFunctionSetVectorFactory > StoredShapeFunctionSetVectorProviderType;

      struct BlockMapperSingletonFactory
      {
        static BlockMapperType *createObject ( LFEMapType *lfeMap )
        {
          return new BlockMapperType( lfeMap->gridPart(), [ lfeMap ] ( const auto &refElement ) {
            if( lfeMap->hasCoefficients( refElement.type() ) )
              return Dune::Fem::compile( refElement, lfeMap->localCoefficients( refElement.type() ) );
            else
              return Dune::Fem::DofMapperCode();
          } );
        }

        static void deleteObject ( BlockMapperType *object ) { delete object; }
      };

      typedef SingletonList< LFEMapType *, BlockMapperType, BlockMapperSingletonFactory > BlockMapperProviderType;

    public:
      typedef LocalFiniteElementInterpolation< ThisType, LocalInterpolationType, Traits::isScalar > InterpolationImplType;
      typedef LocalFEInterpolationWrapper< ThisType > InterpolationType;

      using BaseType::order;

      template< class GridPart, std::enable_if_t< std::is_same< GridPart, GridPartType >::value &&std::is_same< KeyType, std::tuple<> >::value, int > = 0 >
      explicit LocalFiniteElementSpace ( GridPart &gridPart,
                                         const InterfaceType commInterface = InteriorBorder_All_Interface,
                                         const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection ),
          lfeMap_( &LFEMapProviderType::getObject( std::make_pair( &gridPart, KeyType() ) ) ),
          storedShapeFunctionSetVector_( &StoredShapeFunctionSetVectorProviderType::getObject( lfeMap_.get() ) ),
          blockMapper_( &BlockMapperProviderType::getObject( lfeMap_.get() ) )
      {}

      template< class GridPart, std::enable_if_t< std::is_same< GridPart, GridPartType >::value && !std::is_same< KeyType, std::tuple<> >::value, int > = 0 >
      explicit LocalFiniteElementSpace ( GridPart &gridPart, const KeyType &key,
                                         const InterfaceType commInterface = InteriorBorder_All_Interface,
                                         const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection ),
          lfeMap_( &LFEMapProviderType::getObject( std::make_pair( &gridPart, key ) ) ),
          storedShapeFunctionSetVector_( &StoredShapeFunctionSetVectorProviderType::getObject( lfeMap_.get() ) ),
          blockMapper_( &BlockMapperProviderType::getObject( lfeMap_.get() ) )
      {}

      template< class GridPart, std::enable_if_t< std::is_same< GridPart, GridPartType >::value && !std::is_same< KeyType, std::tuple<> >::value, int > = 0 >
      explicit LocalFiniteElementSpace ( GridPart &gridPart,
                                         const InterfaceType commInterface = InteriorBorder_All_Interface,
                                         const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection ),
          lfeMap_( &LFEMapProviderType::getObject( std::make_pair( &gridPart, KeyType(1) ) ) ),
          storedShapeFunctionSetVector_( &StoredShapeFunctionSetVectorProviderType::getObject( lfeMap_.get() ) ),
          blockMapper_( &BlockMapperProviderType::getObject( lfeMap_.get() ) )
      {}
      LocalFiniteElementSpace ( const ThisType & ) = delete;
      LocalFiniteElementSpace ( ThisType && ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type **/
      DFSpaceIdentifier type () const { return LocalFiniteElementSpace_id; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::basisFunctionSet **/
      BasisFunctionSetType basisFunctionSet ( const EntityType &entity ) const
      {
        return BasisFunctionSetType( entity, shapeFunctionSet( entity ) );
      }

      /**
       * \brief return shape function set for given entity
       *
       * \param[in]  entity  entity for which shape function set is requested
       *
       * \returns  ShapeFunctionSetType  shape function set
       **/
      ShapeFunctionSetType shapeFunctionSet ( const EntityType &entity ) const
      {
        return getShapeFunctionSet( (*lfeMap_)( entity ), entity.type() );
      }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous () const { return true; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous ( const IntersectionType &intersection ) const { return true; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::order */
      int order () const { return lfeMap_->order(); }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::multipleGeometryTypes */
      bool multipleGeometryTypes () const { return true; }

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::blockMapper */
      BlockMapperType &blockMapper () const { assert( blockMapper_ ); return *blockMapper_; }

      /**
       * \brief return local interpolation
       *
       *  \param[in]  entity  grid part entity
       **/
      InterpolationType interpolation () const
      {
        return InterpolationType( *this );
      }

      /**
       * \brief return local interpolation
       *
       *  \param[in]  entity  grid part entity
       **/
      [[deprecated("Use LocalInterpolation( space ) instead!")]]
      InterpolationImplType interpolation ( const EntityType &entity ) const
      {
        return localInterpolation( entity );
      }

      /**
       * \brief return local interpolation
       *
       *  \param[in]  entity  grid part entity
       **/
      InterpolationImplType localInterpolation ( const EntityType &entity ) const
      {
        auto lfe = (*lfeMap_)( entity );
        return InterpolationImplType( BasisFunctionSetType( entity, getShapeFunctionSet( lfe, entity.type() ) ), std::get< 2 >( lfe ) );
      }

      typedef typename LFEMapType::LocalCoefficientsType QuadratureType;
      const QuadratureType& quadrature ( const GeometryType &type ) const
      {
        return (*lfeMap_).localCoefficients(type);
      }

    private:
      ShapeFunctionSetType getShapeFunctionSet ( std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & > lfe, const GeometryType &type ) const
      {
        auto &storedShapeFunctionSet = (*storedShapeFunctionSetVector_)[ std::get< 0 > ( lfe ) ];
        if( !storedShapeFunctionSet )
          storedShapeFunctionSet.reset( new StoredShapeFunctionSetType( type, LocalFunctionsShapeFunctionSetType( std::get< 1 >( lfe ) ) ) );
        return ShapeFunctionSetType( storedShapeFunctionSet.get() );
      }

      std::unique_ptr< LFEMapType, typename LFEMapProviderType::Deleter > lfeMap_;
      std::unique_ptr< StoredShapeFunctionSetVectorType, typename StoredShapeFunctionSetVectorProviderType::Deleter > storedShapeFunctionSetVector_;
      std::unique_ptr< BlockMapperType, typename BlockMapperProviderType::Deleter > blockMapper_;
    };

    template< class LFEMap, class FunctionSpace, class Storage, int newRange >
    struct ToNewDimRangeFunctionSpace<
      LocalFiniteElementSpace<LFEMap, FunctionSpace, Storage>, newRange>
    {
      typedef LocalFiniteElementSpace<LFEMap, typename ToNewDimRangeFunctionSpace<FunctionSpace,newRange>::Type, Storage> Type;
    };
    template <class LFEMap, class FunctionSpace, class Storage, class NewFunctionSpace>
    struct DifferentDiscreteFunctionSpace<
        LocalFiniteElementSpace<LFEMap,FunctionSpace,Storage>, NewFunctionSpace>
    {
      typedef LocalFiniteElementSpace<LFEMap, NewFunctionSpace, Storage > Type;
    };
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_SPACE_HH
