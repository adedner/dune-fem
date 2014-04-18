#ifndef DUNE_FEM_TEST_TESTGRID_HH
#define DUNE_FEM_TEST_TESTGRID_HH

#include <sstream>

namespace Dune
{

  class TestGrid
  {
  private:
    typedef TestGrid ThisType;
    typedef GridSelector::GridType GridType;
  protected:
    GridPtr< GridType > gridptr_;
    
  protected:
    inline TestGrid ()
    : gridptr_( macroGridName() )
    {
    }

  private:
    TestGrid ( const ThisType & );

    ThisType &operator= ( const ThisType & );

  public:
    static inline ThisType &instance ()
    {
      static ThisType staticInstance;
      return staticInstance;
    }

    static inline GridType &grid ()
    {
      return *(instance().gridptr_);
    }

    static inline int refineStepsForHalf ()
    {
      return DGFGridInfo< GridType > :: refineStepsForHalf();
    }

  protected:
    static inline std :: string macroGridName ()
    {
      std :: ostringstream s;
      s << GridType :: dimension << "dgrid.dgf";
      return s.str();
    }
  };
  
}

#endif