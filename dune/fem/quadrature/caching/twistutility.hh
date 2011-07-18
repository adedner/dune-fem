#ifndef DUNE_TWISTUTILITY_HH
#define DUNE_TWISTUTILITY_HH

#include <cassert>

#include <dune/common/static_assert.hh>
#include <dune/common/version.hh>
#include <dune/common/geometrytype.hh>
#include <dune/grid/genericgeometry/topologytypes.hh>

#include <dune/grid/alugrid/common/interfaces.hh>

#if HAVE_DUNE_GEOGRID
#include <dune/grid/utility/hostgridaccess.hh>
#endif

#include <dune/fem/misc/capabilities.hh>

namespace Dune
{

  /** \brief Utility to get twist from IntersectionIterator, 
      if provided by grid (i.e. AlbertaGrid, ALUGrid)
      otherwise return default values (correct for YASP/SGRID).
      
      The twist (t) of a face is defined in the following way:
      - sign(t) gives information on the relationship between the
        orientation of the intersection geometry and the geometry
        of the corresponding codim 1 entity of the inside/outside
        entity:
        - sign(t)>=0: same orientation
        - sign(t)<0:  opposite orientation

      - The value of the twist gives information on the local numbering
        of the corners of the corresponding geometries. This value
        is only correctly defined for conforming grids, i.e.,
        the intersection is identical to an codim 1 entity of inside/outside.
        In this case we have the following definition:
        - sign(t)>=0: corner[0] of inside/outside face is equal to
                      corner[t] of intersection.
        - sign(t)<0:  corner[0] of inside/outside face is equal to
                      corner[t'] of intersection with t' = abs(t)+1.
  */
  template< class Grid > 
  class TwistUtility
  {
#if 0
    dune_static_assert( (!Conversion< Grid, HasHierarchicIndexSet > :: exists),
                        "The default TwistUtility is only for SGrid, "
			"YaspGrid, UGGrid and OneDGrid." );
#endif

  public:
    typedef Grid GridType;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) {}

    //! \brief return 0 for inner face 
    template <class IntersectionIterator> 
    static inline int twistInSelf(const GridType &, const IntersectionIterator&)
    {
      return 0;
    }
    
    //! \brief return 0 for inner face 
    template <class IntersectionIterator> 
    int twistInSelf(const IntersectionIterator& it) const {
      return 0;
    }
    
    //! \brief return 0 for outer face 
    template <class IntersectionIterator> 
    int twistInNeighbor(const IntersectionIterator& it) const {
      return 0;
    }

    //! \brief return 0 for outer face 
    template <class IntersectionIterator> 
    static inline int twistInNeighbor(const GridType &, const IntersectionIterator&) 
    {
      return 0;
    }

    /** \brief return geometry type of inside or outside entity */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside) 
    {
      if( Capabilities::IsUnstructured<GridType>::v ) 
      {
        return (inside) ? intersection.inside()->type() :  
                          intersection.outside()->type();  
      }
      else 
      {
#ifndef NDEBUG
        GeometryType geoType( GenericGeometry::CubeTopology< GridType :: dimension >::type::id, 
                              GridType :: dimension );
        GeometryType realType = (inside) ? intersection.inside()->type() :
                                           intersection.outside()->type();
        assert ( realType == geoType );
#endif
        return GeometryType( GenericGeometry::CubeTopology< GridType :: dimension >::type::id, 
                             GridType :: dimension );
      }
    }
  };



  // Specialization for AlbertaGrid
  // ------------------------------
  
#ifdef ENABLE_ALBERTA
  template< int dim, int dimW >
  class AlbertaGrid;
  
  /** \brief Specialization of TwistUtility for AlbertaGrid. 
  */
  template< int dim, int dimW >
  struct TwistUtility< AlbertaGrid< dim, dimW > >
  {
    typedef AlbertaGrid<dim, dimW> GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator  LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection  LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    static const int dimension = GridType::dimension;

  private:
    const GridType &grid_;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility ( const GridType &grid )
    : grid_( grid )
    {}

    //! \brief return twist for inner face 
    static int twistInSelf ( const GridType &grid, const LeafIntersection &it )
    {
      return grid.getTwistInInside( it );
    }
    
    //! \brief return twist for inner face 
    int twistInSelf ( const LeafIntersection &it ) const
    {
      return twistInSelf( grid_, it );
    }
    
    //! \brief return twist for outer face 
    static int twistInNeighbor ( const GridType &grid, const LeafIntersection &it )
    {
      return grid.getTwistInOutside( it );
    }

    //! \brief return twist for outer face 
    int twistInNeighbor ( const LeafIntersection &it ) const
    {
      return twistInNeighbor( grid_, it );
    }
    
    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GenericGeometry::SimplexTopology< dimension >::type::id,
                           dimension );
    }
  };
#endif



  // Specialization for ALUGrid
  // --------------------------

#ifdef ENABLE_ALUGRID
  template< int dim, int dimW >
  class ALUSimplexGrid;

  template< int dim, int dimW >
  class ALUCubeGrid;

  template< int dim, int dimW >
  class ALUConformGrid;

  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template< int dim >
  struct TwistUtility< ALUSimplexGrid< dim, dim > >
  {
    typedef ALUSimplexGrid< dim, dim > GridType;

    typedef typename GridType::Traits::LeafIntersectionIterator
      LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;

    typedef typename GridType::Traits::LevelIntersectionIterator
      LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    static const int dimension = GridType::dimension;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    //! \brief return twist for inner face 
    template <class Intersection> 
    static inline int twistInSelf(const GridType & grid, 
                                  const Intersection& intersection)
    {
      return grid.getRealIntersection( intersection ).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersection& intersection) const 
    {
      return twistInSelf( grid_, intersection );
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersection& intersection) const {
      return twistInSelf( grid_, intersection );
    }

    //! \brief return twist for outer face 
    template <class Intersection>
    static inline int twistInNeighbor(const GridType & grid, 
                                      const Intersection& intersection) 
    {
      return grid.getRealIntersection( intersection ).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersection& intersection) const 
    {
      return twistInNeighbor( grid_, intersection );
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersection& it) const 
    {
      return twistInNeighbor( grid_, it );
    }

    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GenericGeometry::SimplexTopology< dim >::type::id,
                           dim );
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
    
  private:
    const GridType& grid_; 
  };

  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template< int dim >
  struct TwistUtility< ALUCubeGrid< dim, dim > >
  {
    typedef ALUCubeGrid< dim, dim > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    //! \brief return twist for inner face 
    template <class Intersection> 
    static inline int twistInSelf(const GridType & grid, 
                                  const Intersection& intersection)
    {
      return grid.getRealIntersection( intersection ).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersection& intersection) const {
      return twistInSelf( grid_, intersection );
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersection& intersection) const {
      return twistInSelf( grid_, intersection );
    }

    //! \brief return twist for outer face 
    template <class Intersection>
    static inline int twistInNeighbor(const GridType & grid, 
                                      const Intersection& intersection) 
    {
      return grid.getRealIntersection( intersection ).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersection& intersection) const 
    {
      return twistInNeighbor( grid_, intersection);
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersection& intersection) const 
    {
      return twistInNeighbor( grid_, intersection);
    }

    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GenericGeometry::CubeTopology< dim >::type::id,
                           dim );
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
    
  private:
    const GridType& grid_; 
  };

  template< int dim >
  struct TwistUtility< ALUConformGrid< dim, dim > >
  {
    typedef ALUConformGrid< dim, dim > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) : grid_(grid) {}

    //! \brief return twist for inner face 
    static inline int twistInSelf(const GridType & grid, const LeafIntersection& intersection)
    {
      return grid.getRealIntersection( intersection ).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersection& intersection) const {
      return twistInSelf( grid_, intersection );
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersection& intersection) const {
      return twistInSelf( grid_, intersection );
    }

    //! \brief return twist for outer face 
    static inline int twistInNeighbor(const GridType &grid, const LeafIntersection& intersection )
    {
      return grid.getRealIntersection( intersection ).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersection& intersection) const {
      return twistInNeighbor( grid_, intersection );
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersection& intersection) const {
      return twistInNeighbor( grid_, intersection );
    }

    /** \brief return element geometry type of inside or outside entity 
    */
    template <class Intersection>  
    static inline GeometryType
    elementGeometry(const Intersection& intersection, 
                    const bool inside)
    {
      return GeometryType( GenericGeometry::SimplexTopology< dim >::type::id,
                           dim );
    }
    
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);

  private:
    const GridType& grid_; 
  };
#endif



  // Specialization for UGGrid
  // -------------------------

#ifdef ENABLE_UG
  template< int dim >
  class UGGrid;

  template< int dim >
  struct TwistUtility< UGGrid< dim > >
  {
    typedef UGGrid< dim > GridType;

    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    explicit TwistUtility ( const GridType &grid )
    {}

    static int twistInSelf ( const GridType &grid, const LeafIntersection &it );
    int twistInSelf( const LeafIntersection &it ) const;
    static int twistInSelf ( const GridType &grid, const LevelIntersection &it );
    int twistInSelf( const LevelIntersection &it ) const;

    static int twistInNeighbor ( const GridType &grid, const LeafIntersection &it );
    int twistInNeighbor( const LeafIntersection &it ) const;
    static int twistInNeighbor ( const GridType &grid, const LevelIntersection &it );
    int twistInNeighbor( const LevelIntersection &it ) const;

    template< class Intersection >
    static GeometryType
    elementGeometry ( const Intersection &intersection, const bool inside )
    {
      return (inside ? intersection.inside()->type() : intersection.outside()->type());
    }
  };
#endif // #ifdef ENABLE_UG



  // Specialization for GeoGrid
  // --------------------------

#if HAVE_DUNE_GEOGRID
  template< class HostGrid, class CoordFunction >
  class GeometryGrid;

  template< class HostGrid, class CoordFunction >
  struct TwistUtility< GeometryGrid< HostGrid, CoordFunction > >
  {
    typedef GeometryGrid< HostGrid, CoordFunction > GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

  private:
    typedef TwistUtility< HostGrid > HostTwistUtility;
    typedef Dune::HostGridAccess< GridType > HostGridAccess;

    const GridType &grid_;

  public:
    //! \brief constructor
    TwistUtility ( const GridType &grid )
    : grid_( grid )
    {}

    //! \brief return twist for inner face
    static int twistInSelf ( const GridType &grid,
                             const LeafIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLeafIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: twistInSelf( grid.hostGrid(), hostIntersection );
    }

    //! \brief return twist for inner face
    static int twistInSelf ( const GridType &grid,
                             const LevelIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLevelIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: twistInSelf( grid.hostGrid(), hostIntersection );
    }

    //! \brief return twist for inner face
    template< class Intersection >
    int twistInSelf ( const Intersection &intersection ) const
    {
      return twistInSelf( grid_, intersection );
    }


    //! \brief return twist for outer face
    static int twistInNeighbor ( const GridType &grid,
                                 const LeafIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLeafIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: twistInNeighbor( grid.hostGrid(), hostIntersection );
    }

    //! \brief return twist for outer face
    static int twistInNeighbor ( const GridType &grid,
                                 const LevelIntersection &intersection )
    {
      typedef typename HostGridAccess :: HostLevelIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: twistInNeighbor( grid.hostGrid(), hostIntersection );
    }

    //! \brief return twist for outer face
    template< class Intersection >
    int twistInNeighbor ( const Intersection &intersection ) const
    {
      return twistInNeighbor( grid_, intersection );
    }

    /** \brief return element geometry type of inside or outside entity */
    static GeometryType
    elementGeometry ( const LeafIntersection &intersection, bool inside )
    {
      typedef typename HostGridAccess :: HostLeafIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: elementGeometry( hostIntersection, inside );
    }

    /** \brief return element geometry type of inside or outside entity */
    static GeometryType
    elementGeometry ( const LevelIntersection &intersection, bool inside )
    {
      typedef typename HostGridAccess :: HostLevelIntersection HostIntersection;
      const HostIntersection &hostIntersection
        = HostGridAccess :: getIntersection( intersection );
      return HostTwistUtility :: elementGeometry( hostIntersection, inside );
    }

  private:
    TwistUtility( const TwistUtility & );
    TwistUtility &operator=( const TwistUtility & );
  };
#endif // #if HAVE_DUNE_GEOGRID
  
} // end namespace Dune 

#endif // #ifndef DUNE_TWISTUTILITY_HH
