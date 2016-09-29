#ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
#define DUNE_FEMPY_PY_GRID_GRIDPART_HH

#include <cassert>

#include <map>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>

#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/grid/gridpartadapter.hh>
#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/py/grid/function.hh>
#include <dune/fempy/py/grid/numpy.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // GridPartConverter
      // -----------------

      template< class GV >
      struct GridPartConverter
      {
        typedef GV GridView;
        typedef GridPartAdapter< GV > GridPart;

        GridPart &operator() ( pybind11::handle gridView )
        {
          auto result = instances_.emplace( gridView.ptr(), nullptr );
          auto pos = result.first;
          if( result.second )
          {
            // create new gridpart object
            pos->second = new GridPart( gridView.template cast< GridView >() );

            // create Python guard object, removing the grid part once the grid view dies
            pybind11::cpp_function remove_gridpart( [ this, pos ] ( pybind11::handle weakref ) {
                delete pos->second;
                instances_.erase( pos );
                weakref.dec_ref();
              } );
            pybind11::weakref weakref( gridView, remove_gridpart );
            weakref.release();
          }
          assert( pos->second );
          return *pos->second;
        }

      private:
        std::map< PyObject *, GridPart * > instances_;
      };


      template< class GP >
      struct GridPartConverter< Dune::GridView< Fem::GridPart2GridViewTraits< GP > > >
      {
        typedef GP GridPart;
        typedef Dune::GridView< Fem::GridPart2GridViewTraits< GP > > GridView;

        GridPart &operator() ( pybind11::handle gridView )
        {
          return const_cast< GridPart & >( gridView.template cast< GridView >().impl().gridPart() );
        }
      };


      template< class GP >
      struct GridPartConverter< Fem::GridPart2GridView< GP > >
        : public GridPartConverter< Dune::GridView< Fem::GridPart2GridViewTraits< GP > > >
      {};



      // gridPartConverter
      // -----------------

      template< class GridView >
      inline GridPartConverter< GridView > &gridPartConverter ()
      {
        static GridPartConverter< GridView > converter;
        return converter;
      }

    } // namespace detail



    // GridPart
    // --------

    template< class GridView >
    using GridPart = typename detail::GridPartConverter< GridView >::GridPart;



    // gridPart
    // --------

    template< class GridView >
    inline static GridPart< GridView > &gridPart ( pybind11::handle gridView )
    {
      return detail::gridPartConverter< GridView >()( std::move( gridView ) );
    }


#if 0
    // registerGridPart
    // ----------------

    namespace detail
    {

      template< class GridPart, class Field, int dimR >
      double l2Norm ( pybind11::object &o )
      {
        typedef Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< Field, dimR > > VF;
        const VF &u = o.cast< VF >();
        Dune::Fem::L2Norm< GridPart > norm( u.gridPart(), 8 );
        return norm.norm( u );
      }

      template< class GridPart, int ... dimRange >
      auto defL2Norm( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange ... > )
      {
        typedef std::function< double (pybind11::object &) > Dispatch;
        std::array< Dispatch, sizeof ... ( dimRange ) > dispatch = {{ Dispatch( l2Norm< GridPart, double, dimRange > ) ... }};

        return [ dispatch ] ( pybind11::object gp, pybind11::object func ) {
           //const GridPart &gridPart = gp.cast< const GridPart & >();
           pybind11::object dimRobj = func.attr( "dimRange" );
           int dimR = func.attr( "dimRange" ).cast< int >();
           if( static_cast< std::size_t >( dimR ) >= dispatch.size() )
             DUNE_THROW( NotImplemented, "localGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
           return dispatch[ static_cast< std::size_t >( dimR ) ]( func );
        };
      }

      template< class GridPart, class Cls >
      void registerGridPartConstructorFromGrid ( Cls &cls, std::false_type )
      {}

      template< class GridPart, class Cls >
      void registerGridPartConstructorFromGrid ( Cls &cls, std::true_type )
      {
        typedef typename GridPart::GridType Grid;
        cls.def( "__init__", [] ( GridPart &instance, Grid &grid ) {
            new (&instance)GridPart( grid );
          }, pybind11::keep_alive< 1, 2 >() );
      }

      template< class GridPart, class Cls >
      void registerGridPart ( pybind11::handle scope, Cls &cls )
      {
        typedef typename GridPart::GridType Grid;

        registerGridPartConstructorFromGrid< GridPart, Cls >( cls, std::is_constructible< GridPart, Grid & >() );
        cls.attr( "dimGrid" ) = pybind11::int_( GridPart::dimension );
        cls.attr( "dimWorld" ) = pybind11::int_( GridPart::dimensionworld );

        cls.def( "globalGridFunction", defGlobalGridFunction< GridPart >( cls, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ) );
        cls.def( "localGridFunction", defLocalGridFunction< GridPart >( cls, "LocalGridFunction", std::make_integer_sequence< int, 11 >() ) );

        cls.def( "l2Norm", defL2Norm< GridPart >( cls, "l2Norm", std::make_integer_sequence< int, 11 >() ) );
      }

    } // namespace detail
#endif // #if 0

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
