#ifndef DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH
#define DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH

#if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
#error "Experimental grid extensions required for BoundaryIdProvider. Reconfigure with -DDUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS=TRUE."
#endif // #if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid/declaration.hh>
#endif // #if HAVE_DUNE_SPGRID

#include <dune/fem/misc/griddeclaration.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Grid >
  struct HostGridAccess;



  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Grid >
    struct BoundaryIdProvider;



    // BoundaryIdProvider for AlbertaGrid
    // ----------------------------------

#if HAVE_ALBERTA
    template< int dim, int dimW >
    struct BoundaryIdProvider< AlbertaGrid< dim, dimW > >
    {
      typedef AlbertaGrid< dim, dimW > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.impl().boundaryId();
      }
    };
#endif // #if HAVE_ALBERTA



    // BoundaryIdProvider for ALUGrid
    // ------------------------------

#if HAVE_DUNE_ALUGRID
    template< int dim, int dimw, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct BoundaryIdProvider< ALUGrid< dim, dimw, elType, refineType, Comm > >
    {
      typedef ALUGrid< dim, dimw, elType, refineType, Comm > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.impl().boundaryId();
      }
    };
#endif // #if HAVE_DUNE_ALUGRID



    // BoundaryIdProvider for CacheItGrid
    // ----------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< CacheItGrid< HostGrid > >
    {
      typedef CacheItGrid< HostGrid > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for CartesianGrid
    // ------------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< CartesianGrid< HostGrid > >
    {
      typedef CartesianGrid< HostGrid > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for FilteredGrid
    // -----------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< FilteredGrid< HostGrid > >
    {
      typedef FilteredGrid< HostGrid > GridType;

      // todo: FilteredGrid is a filtering grid and, hence, needs a specialized
      //       version of boundaryId.
      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        if( !HostGridAccess< GridType >::hostIntersection( intersection ).boundary() )
          DUNE_THROW( NotImplemented, "BoundaryIdProvider for artificial boundaries of FilteredGrid not implemented." );
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for GeometryGrid
    // -----------------------------------

    template< class HostGrid, class CoordFunction, class Allocator >
    struct BoundaryIdProvider< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      typedef GeometryGrid< HostGrid, CoordFunction, Allocator > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };



    // BoundaryIdProvider for IdGrid
    // -----------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< IdGrid< HostGrid > >
    {
      typedef IdGrid< HostGrid > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for OneDGrid
    // -------------------------------

    template<>
    struct BoundaryIdProvider< OneDGrid >
    {
      typedef OneDGrid GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.boundarySegmentIndex();
      }
    };



    // BoundaryIdProvider for ParallelGrid
    // -----------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct BoundaryIdProvider< ParallelGrid< HostGrid > >
    {
      typedef ParallelGrid< HostGrid > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for SphereGrid
    // ---------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid, class MapToSphere >
    struct BoundaryIdProvider< SphereGrid< HostGrid, MapToSphere > >
    {
      typedef SphereGrid< HostGrid, MapToSphere > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return BoundaryIdProvider< HostGrid >
          ::boundaryId ( HostGridAccess< GridType >::hostIntersection( intersection ) );
      }
    };
#endif // #if HAVE_DUNE_METAGRID



    // BoundaryIdProvider for SPGrid
    // -----------------------------

#if HAVE_DUNE_SPGRID
    template< class ct, int dim, template< int > class Strategy, class Comm >
    struct BoundaryIdProvider< SPGrid< ct, dim, Strategy, Comm > >
    {
      typedef SPGrid< ct, dim, Strategy, Comm > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.indexInInside()+1) : 0);
      }
    };
#endif // #if HAVE_DUNE_SPGRID



    // BoundaryIdProvider for UGGrid
    // -----------------------------

    template< int dim >
    struct BoundaryIdProvider< UGGrid< dim > >
    {
      typedef UGGrid< dim > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return intersection.boundarySegmentIndex();
      }
    };



    // BoundaryIdProvider for YaspGrid
    // -------------------------------

    template< int dim, class CoordCont >
    struct BoundaryIdProvider< YaspGrid< dim, CoordCont > >
    {
      typedef YaspGrid< dim, CoordCont > GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.indexInInside()+1) : 0);
      }
    };

    ///////////////////////////////////////////////////////

    template< class Grid >
    struct BoundaryIdGetter
    {
      typedef std::function<int(const typename Grid::template Codim<0>::Geometry::GlobalCoordinate&)> Caller;
      BoundaryIdGetter(const Caller &caller)
      : caller_(caller) {}
      template< class Intersection >
      int boundaryId ( const Intersection &intersection ) const
      {
        if (!intersection.boundary) return 0;
        const auto x = intersection.geometry().center();
        int val = caller(x);
        return (val==0)?BoundaryIdProvider<Grid>::boundaryId(intersection):val;
      }
      private:
      Caller caller_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_BOUNDARYIDPROVIDER_HH
