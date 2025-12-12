#ifndef DUNE_FEM_GRIDPART_DUALGRIDPART_HH
#define DUNE_FEM_GRIDPART_DUALGRIDPART_HH

#include <dune/grid/common/gridview.hh>

#include <dune/fem/gridpart/common/deaditerator.hh>
#include <dune/fem/gridpart/common/entitysearch.hh>
#include <dune/fem/gridpart/common/extendedentity.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/space/common/dofmanager.hh>

//#include <dune/fem/gridpart/dualgridpart/capabilities.hh>

#include <dune/polygongrid/grid.hh>
//#include <dune/fem/gridpart/dualgridpart/datahandle.hh>
//#include <dune/fem/gridpart/dualgridpart/entity.hh>
//#include <dune/fem/gridpart/dualgridpart/geometry.hh>
//#include <dune/fem/gridpart/dualgridpart/indexset.hh>
//#include <dune/fem/gridpart/dualgridpart/intersection.hh>
//#include <dune/fem/gridpart/dualgridpart/intersectioniterator.hh>
//#include <dune/fem/gridpart/dualgridpart/iterator.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class HostGridPartImp >
    class DualGridPart;


    // DualGridPartTraits
    // ------------------

    template< class HostGridPartImp >
    struct DualGridPartTraits
    {
      typedef DualGridPart< HostGridPartImp > GridPartType;

      typedef HostGridPartImp  HostGridPartType;

      typedef typename HostGridPartType::ctype ctype;
      typedef Dune::PolygonGrid< ctype > GridType;

      /** \brief The type of the corresponding TwistUtility */
      typedef TwistUtility< GridType >  TwistUtilityType ;
      //! type of twist utility
      //typedef MetaTwistUtility< PolyTwistUtilityType >  TwistUtilityType;

      typedef typename GridType::LeafGridView LeafGridView;

      //! default partition iterator type
      static const PartitionIteratorType indexSetPartitionType = HostGridPartType::indexSetPartitionType;
      static const InterfaceType indexSetInterfaceType = HostGridPartType::indexSetInterfaceType;


      template< int codim >
      struct Codim
      {
        typedef typename GridType::Traits::template Codim<codim>::Entity EntityType;
        typedef EntityType Entity;

        typedef typename GridType::Traits::template Codim<codim>::Geometry GeometryType;
        typedef GeometryType Geometry;

        typedef typename GridType::Traits::template Codim<codim>::LocalGeometry LocalGeometryType;
        typedef LocalGeometryType LocalGeometry;

        typedef typename GridType::Traits::template Codim<codim>::EntitySeed   EntitySeed;
        typedef EntitySeed  EntitySeedType;


        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef typename LeafGridView::template Codim<codim>::template Partition<pitype>::Iterator Iterator;
          typedef Iterator  IteratorType;
        };
      };

      typedef typename LeafGridView::IndexSet  IndexSetType;

      typedef typename LeafGridView::IntersectionIterator  IntersectionIterator;
      typedef IntersectionIterator IntersectionIteratorType;

      typedef typename GridType::Communication CommunicationType;

      static const bool conforming = true;
    };



    // DualGridPart
    // ------------

    template< class HostGridPartImp >
    class DualGridPart
    : public GridPartDefault< DualGridPartTraits< HostGridPartImp > >
    {
      typedef DualGridPart< HostGridPartImp > ThisType;
      typedef GridPartDefault< DualGridPartTraits< HostGridPartImp > > BaseType;


    public:
      typedef DualGridPartTraits< HostGridPartImp > Traits;

      typedef typename Traits::HostGridPartType HostGridPartType;

      typedef typename BaseType::GridType GridType;

      typedef typename Traits::IndexSetType IndexSetType;
      typedef typename BaseType::IntersectionIteratorType IntersectionIteratorType;
      typedef typename BaseType::IntersectionType IntersectionType;
      typedef typename BaseType::CommunicationType CommunicationType;

      using BaseType::grid;

      // from dune-polygongrid
      typedef typename GridType::Mesh      Mesh;
      typedef typename GridType::MeshType  MeshType;


      typedef DofManager< GridType > DofManagerType;
      typedef DofManager< typename HostGridPartType::GridType > HostDofManagerType;

      template< int codim >
      struct Codim
      : public BaseType::template Codim< codim >
      {};

      explicit DualGridPart ( typename HostGridPartType::GridType &og )
      : BaseType( *createGrid( HostGridPartType( og ) ) ),
        gridPtr_( &this->grid() ),
        hostGridPart_( og ),
        leafGridView_( gridPtr_->leafGridView() ),
        dofManager_( &DofManagerType :: instance( this->grid() ) ),
        hostDofManager_( &HostDofManagerType :: instance( hostGridPart_.grid() ) ),
        sequence_( -1 )
      {}

      explicit DualGridPart ( const DualGridPart &other )
      : BaseType( other ),
        gridPtr_( other.gridPtr_ ),
        hostGridPart_( other.hostGridPart() ),
        leafGridView_( gridPtr_->leafGridView() ),
        dofManager_( &DofManagerType :: instance( this->grid() ) ),
        hostDofManager_( &HostDofManagerType :: instance( hostGridPart_.grid() ) ),
        sequence_( -1 )
      {}

      DualGridPart& operator= ( const DualGridPart& other ) = default;

      explicit DualGridPart ( const HostGridPartType &hostGridPart )
      : BaseType( *createGrid( hostGridPart ) ),
        gridPtr_( &this->grid() ),
        hostGridPart_( hostGridPart ),
        leafGridView_( gridPtr_->leafGridView() ),
        dofManager_( &DofManagerType :: instance( this->grid() ) ),
        hostDofManager_( &HostDofManagerType :: instance( hostGridPart_.grid() ) ),
        sequence_( -1 )
      {
      }

      void update() const
      {
        const int hostSequence = hostDofManager_->sequence();
        if( sequence_ != hostSequence )
        {
          gridPtr_->update( createMesh( hostGridPart_ ), __PolygonGrid::Dual );
          dofManager_->resizeMemory();
          sequence_ = hostSequence;
        }
      }

      const IndexSetType &indexSet () const
      {
        update();
        return leafGridView_.indexSet();
      }

      template< int codim >
      typename Codim< codim >::IteratorType
      begin () const
      {
        update();
        return leafGridView_.template begin< codim >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      begin () const
      {
        update();
        return leafGridView_.template begin< codim, pitype >();
      }

      template< int codim >
      typename Codim< codim >::IteratorType
      end () const
      {
        return leafGridView_.template end< codim >();
      }

      template< int codim, PartitionIteratorType pitype >
      typename Codim< codim >::template Partition< pitype >::IteratorType
      end () const
      {
        return leafGridView_.template end< codim, pitype >();
      }

      int level () const
      {
        return 0;
      }

      IntersectionIteratorType ibegin ( const typename Codim< 0 >::EntityType &entity ) const
      {
        update();
        return leafGridView_.ibegin( entity );
      }

      IntersectionIteratorType iend ( const typename Codim< 0 >::EntityType &entity ) const
      {
        update();
        return leafGridView_.iend( entity );
      }

      template< class DataHandle, class Data >
      void communicate ( CommDataHandleIF< DataHandle, Data > &handle,
                         InterfaceType iftype, CommunicationDirection dir ) const
      {
        //typedef CommDataHandleIF< DataHandle, Data >  HostHandleType;
        //IdDataHandle< HostHandleType, GridFamily > handleWrapper( data(), handle );
        //hostGridPart().communicate( handleWrapper, iftype, dir );
      }

      template < class EntitySeed >
      typename Codim< EntitySeed::codimension >::EntityType
      entity ( const EntitySeed &seed ) const
      {
        return grid().entity( seed );
      }

      // convert a grid entity to a grid part entity ("Gurke!")
      template< class Entity >
      typename Codim< Entity::codimension >::EntityType
      //MakeableInterfaceObject< typename Codim< Entity::codimension >::EntityType >
      convert ( const Entity &entity ) const
      {
        return entity;
      }

      const HostGridPartType &hostGridPart () const { return hostGridPart_; }
      HostGridPartType &hostGridPart () { return hostGridPart_; }

    protected:
      std::shared_ptr< Mesh > createMesh(const HostGridPartType& hostGridPart) const
      {
        const int dim = HostGridPartType::dimension;
        const auto& hostIdxSet = hostGridPart.indexSet();
        typedef typename Codim< 0 >::EntityType EntityType;
        typedef typename EntityType::Geometry::GlobalCoordinate   GlobalCoordinate;
        std::vector< GlobalCoordinate > vertices( hostIdxSet.size(dim) );

        std::vector< std::vector< size_t > > polygons( hostIdxSet.size(0) );

        const auto end = hostGridPart.template end<0>();
        for( auto it = hostGridPart.template begin<0>(); it != end; ++it )
        {
          const auto& entity = *it;
          const auto& geom = entity.geometry();
          const int nVx = geom.corners();
          const int idx = hostIdxSet.index( entity );

          auto& elempoly = polygons[ idx ];
          elempoly.resize( nVx );

          // vertices
          for( int i=0; i<nVx; ++i )
          {
            elempoly[ i ] = hostIdxSet.subIndex( entity, i, dim );
            vertices[ elempoly[ i ] ] = geom.corner( i );
          }

          if( nVx == 4 )
            std::swap( elempoly[ 2 ], elempoly[ 3 ] );

          // insert polygon oriented counter-clockwise
          const GlobalCoordinate a = vertices[ elempoly[ 1 ] ] - vertices[ elempoly[ 0 ] ];
          const GlobalCoordinate b = vertices[ elempoly[ 2 ] ] - vertices[ elempoly[ 0 ] ];
          if( a[ 0 ]*b[ 1 ] < a[ 1 ]*b[ 0 ] )
            std::reverse( elempoly.begin(), elempoly.end() );
        }

        const size_t nPoly = hostIdxSet.size(0);
        std::vector< size_t > counts( nPoly );
        for( size_t i=0; i<nPoly; ++i )
        {
          counts[ i ] = polygons[i].size();
        }

        __PolygonGrid::MultiVector< size_t > polyvec( counts );
        for( size_t i=0; i<nPoly; ++i )
        {
          std::copy( polygons[i].begin(), polygons[i].end(), polyvec[ i ].begin() );
        }
        return std::make_shared< Mesh >(vertices, polyvec );
      }

      GridType* createGrid(const HostGridPartType& hostGridPart) const
      {
        return new GridType( createMesh( hostGridPart ), __PolygonGrid::Dual );
      }

      std::shared_ptr< GridType > gridPtr_;
      HostGridPartType hostGridPart_;
      typename GridType::LeafGridView leafGridView_;

      DofManagerType* dofManager_;
      HostDofManagerType* hostDofManager_;
      mutable int sequence_;
    };



    // GridEntityAccess for DualEntity
    // -------------------------------

    template< int codim, int dim, class GridFamily >
    struct GridEntityAccess< Dune::ExtendedEntity< codim, dim, GridFamily, __PolygonGrid::Entity> >
    {
      typedef Dune::ExtendedEntity< codim, dim, GridFamily, __PolygonGrid::Entity > EntityType;
      typedef Dune::Entity< codim, dim, GridFamily, __PolygonGrid::Entity > GridEntityType;
      //typedef GridEntityAccess< typename GridFamily::Traits::HostGridPartType::Traits::template Codim<codim>::EntityType > HostAccessType;
      //typedef typename HostAccessType::GridEntityType GridEntityType;

      static const GridEntityType &gridEntity ( const EntityType &entity )
      {
        return entity;
      }
    };



    // EntitySearch for DualGridPart
    // ---------------------------

#if 0
    template< class HostGridPart, int codim, PartitionIteratorType partition >
    class EntitySearch< DualGridPart< HostGridPart >, codim, partition >
    {
      typedef EntitySearch< DualGridPart< HostGridPart >, codim, partition > ThisType;

    public:
      typedef DualGridPart< HostGridPart > GridPartType;
      typedef typename GridPartType::ExtraData  ExtraData;

      typedef typename GridPartType::template Codim< codim >::EntityType EntityType;

      typedef typename EntityType::Geometry::GlobalCoordinate GlobalCoordinateType;

      explicit EntitySearch ( const GridPartType &gridPart )
      : hostEntitySearch_( gridPart.hostGridPart() ),
        data_( gridPart.data() )
      {}

      EntityType operator() ( const GlobalCoordinateType &x ) const
      {
        typedef typename EntityType::Implementation EntityImpl;
        return EntityImpl( data_, hostEntitySearch_( x ) );
      }

    protected:
      const EntitySearch< HostGridPart > hostEntitySearch_;
      ExtraData data_;
    };
#endif

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_DUALGRIDPART_HH
