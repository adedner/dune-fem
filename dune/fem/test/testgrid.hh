#ifndef DUNE_FEM_TEST_TESTGRID_HH
#define DUNE_FEM_TEST_TESTGRID_HH

// C++ includes
#include <sstream>

// dune-grid includes
#include <dune/grid/io/file/dgfparser/dgfparser.hh>


namespace Dune
{

  namespace Fem
  {

    // TestGrid
    // --------

    class TestGrid
    {
      typedef TestGrid ThisType;
      typedef Dune::GridSelector::GridType HGridType;

    protected:
      TestGrid ()
      : gridptr_( macroGridName() )
      {
        gridptr_->loadBalance();
      }

    private:
      TestGrid ( const ThisType & );

      ThisType &operator= ( const ThisType & );

    public:
      static ThisType &instance ()
      {
        static ThisType staticInstance;
        return staticInstance;
      }

      static HGridType &grid ()
      {
        return *(instance().gridptr_);
      }

      static int refineStepsForHalf ()
      {
        return DGFGridInfo< HGridType >::refineStepsForHalf();
      }

    protected:
      static std::string macroGridName ()
      {
        std::ostringstream s;
        s << HGridType::dimension << "dgrid.dgf";
        return s.str();
      }

      GridPtr< HGridType > gridptr_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_TEST_TESTGRID_HH
