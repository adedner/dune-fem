#ifndef DUNE_TWISTUTILITY_HH
#define DUNE_TWISTUTILITY_HH

// is Alberta was found then also include headers 
#ifndef HAVE_ALBERTA
#define HAVE_ALBERTA_FOUND 0 
#else 
#define HAVE_ALBERTA_FOUND HAVE_ALBERTA 
#endif

// is ALU3dGrid was found then also include headers 
#ifndef HAVE_ALUGRID
#define HAVE_ALUGRID_FOUND 0 
#else 
#define HAVE_ALUGRID_FOUND HAVE_ALUGRID 
#endif

#if HAVE_ALUGRID_FOUND
#include <dune/grid/alugrid.hh>
#endif

#if HAVE_ALBERTA_FOUND
#include <dune/grid/albertagrid.hh>
#endif

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

namespace Dune {

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
  template <class GridImp> 
  class TwistUtility
  {
    // this default implementation only is for SGrid, YaspGrid, UGGrid
    // and OneDGrid. 
    CompileTimeChecker<Conversion<GridImp,HasHierarchicIndexSet>::exists == false>
      implement_specialized_twist_utility;
  public:
    typedef GridImp GridType;
  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) {}

    //! \brief return 0 for inner face 
    template <class IntersectionIterator> 
    static int twistInSelf(const GridType &, const IntersectionIterator&)
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
    static int twistInNeighbor(const GridType &, const IntersectionIterator&) 
    {
      return 0;
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const { return true; }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool conforming (const GridType &, const IntersectionIterator&) { return true; }
  };
  
#if HAVE_ALBERTA_FOUND
  /** \brief Specialization of TwistUtility for AlbertaGrid. 
  */
  template <int dim, int dimW>
  class TwistUtility<AlbertaGrid<dim, dimW> >
  {
  public:
    typedef AlbertaGrid<dim, dimW> GridType;
    typedef typename GridType::Traits::LeafIntersectionIterator  LeafIntersectionIterator;
    typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    //! \brief return twist for inner face 
    static int twistInSelf(const GridType & grid, 
                           const LeafIntersectionIterator& it) {
      return grid.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //int twistInSelf(const LevelIntersectionIterator& it) const {
    //  return grid_.getRealIntersectionIterator(it).twistInSelf();
    //}

    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    static int twistInNeighbor(const GridType & grid, 
                               const LeafIntersectionIterator& it)
    {
      return grid.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //int twistInNeighbor(const LevelIntersectionIterator& it) const {
    //  return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    //}

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const { return true; }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool conforming (const GridType & grid, 
                            const IntersectionIterator& it) { return true; }
  private:
    const GridType& grid_;
  };
#endif

#if HAVE_ALUGRID_FOUND
  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template <>
  class TwistUtility<ALUSimplexGrid<3,3>  >
  {
  public:
    typedef ALUSimplexGrid<3,3> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    //! \brief return twist for inner face 
    template <class IntersectionIterator> 
    static int twistInSelf(const GridType & grid, 
                           const IntersectionIterator& it)
    {
      return grid.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

    //! \brief return twist for outer face 
    template <class IntersectionIterator>
    static int twistInNeighbor(const GridType & grid, 
                               const IntersectionIterator& it) 
    {
      return grid.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const 
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool conforming (const GridType & grid, 
                     const IntersectionIterator& it)
    { 
      return grid.getRealIntersectionIterator(it).conforming(); 
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
    
  private:
    const GridType& grid_; 
  };

  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template <>
  class TwistUtility<ALUCubeGrid<3,3>  >
  {
  public:
    typedef ALUCubeGrid<3,3> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) :
      grid_(grid)
    {}

    //! \brief return twist for inner face 
    template <class IntersectionIterator> 
    static int twistInSelf(const GridType & grid, 
                           const IntersectionIterator& it)
    {
      return grid.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInSelf();
    }

    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    template <class IntersectionIterator>
    static int twistInNeighbor(const GridType & grid, 
                               const IntersectionIterator& it) 
    {
      return grid.getRealIntersectionIterator(it).twistInNeighbor();
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersectionIterator& it) const {
      return grid_.getRealIntersectionIterator(it).twistInNeighbor();
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const 
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool conforming (const GridType & grid, 
                            const IntersectionIterator& it)
    { 
      return grid.getRealIntersectionIterator(it).conforming(); 
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);
    
  private:
    const GridType& grid_; 
  };
  
  /** \brief Specialization of TwistUtility for ALUGridSimplex. 
  */
  template <>
  class TwistUtility<ALUSimplexGrid<2,2>  >
  {
  public:
    typedef ALUSimplexGrid<2, 2> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) : grid_(grid) {}

    //! \brief return twist for inner face 
    static int twistInSelf(const GridType &, const LeafIntersectionIterator&)
    {
      return 0;
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersectionIterator& it) const {
      return 0;
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersectionIterator& it) const {
      return 0;
    }

    //! \brief return twist for outer face 
    static int twistInNeighbor(const GridType &, const LeafIntersectionIterator&)
    {
      return 1;
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return 1;
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersectionIterator& it) const {
      return 1;
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const 
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool conforming (const GridType & grid, 
                     const IntersectionIterator& it)
    { 
      return grid.getRealIntersectionIterator(it).conforming(); 
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);

  private:
    const GridType& grid_; 
  };
  
  template <>
  class TwistUtility<ALUConformGrid<2,2>  >
  {
  public:
    typedef ALUConformGrid<2, 2> GridType;
    typedef GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
    typedef GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
  public:
    //! \brief constructor taking grid reference 
    TwistUtility(const GridType& grid) : grid_(grid) {}

    //! \brief return twist for inner face 
    static int twistInSelf(const GridType &, const LeafIntersectionIterator&)
    {
      return 0;
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LeafIntersectionIterator& it) const {
      return 0;
    }
    
    //! \brief return twist for inner face 
    int twistInSelf(const LevelIntersectionIterator& it) const {
      return 0;
    }

    //! \brief return twist for outer face 
    static int twistInNeighbor(const GridType &, const LeafIntersectionIterator&)
    {
      return 1;
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LeafIntersectionIterator& it) const {
      return 1;
    }
    
    //! \brief return twist for outer face 
    int twistInNeighbor(const LevelIntersectionIterator& it) const {
      return 1;
    }

    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    bool conforming (const IntersectionIterator& it) const 
    { 
      return grid_.getRealIntersectionIterator(it).conforming(); 
    }
    
    //! \brief return true if intersection is conform, default is true  
    template <class IntersectionIterator> 
    static bool conforming (const GridType & grid, 
                     const IntersectionIterator& it)
    { 
      return grid.getRealIntersectionIterator(it).conforming(); 
    }
  private:
    TwistUtility(const TwistUtility&);
    TwistUtility& operator=(const TwistUtility&);

  private:
    const GridType& grid_; 
  };
#endif

#undef HAVE_ALBERTA_FOUND
#undef HAVE_ALUGRID_FOUND
} // end namespace Dune 

#endif
