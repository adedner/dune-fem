#undef NDEBUG

#include <config.h>
#include <iostream>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/misc/gridwidth.hh>

#if defined USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#elif defined USE_VECTORFUNCTION
#include <dune/fem/storage/vector.hh>
#include <dune/fem/function/vectorfunction.hh>
#elif defined USE_ATTACHEDFUNCTION
#include <dune/fem/function/attachedfunction/function.hh>
#elif defined USE_BLOCKVECTORDISCRETEFUNCTION
#include <dune/fem/function/blockvectordiscretefunction.hh>
#include <dune/fem/function/blockvectors/referenceblockvector.hh>
#elif defined USE_COMBINEDFUNCTION
#undef DIMRANGE 
#define DIMRANGE 1
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/combinedfunction.hh>
#else
#include <dune/fem/function/adaptivefunction.hh>
#endif

#if defined  USE_FILTEREDGRID 
#include <dune/fem/gridpart/filteredgrid.hh>
#endif
#if defined  USE_IDGRIDPART 
#include <dune/fem/gridpart/idgridpart.hh>
#endif

#include "testgrid.hh"
#include "dfspace.hh"
#include "exactsolution.hh"

using namespace Dune;

typedef GridSelector::GridType MyGridType;
// typedef HierarchicGridPart< MyGridType >  ContainedGridPartType;
typedef DGAdaptiveLeafGridPart< MyGridType > ContainedGridPartType;
//typedef IdBasedLeafGridPart< MyGridType > ContainedGridPartType;

// use filtered grid for testing 
#if defined  USE_FILTEREDGRID 
  typedef RadialFilter< ContainedGridPartType > FilterType;
  typedef FilteredGridPart<ContainedGridPartType, FilterType, true > GridPartType;
#elif defined  USE_IDGRIDPART
  typedef Fem:: IdGridPart< ContainedGridPartType > GridPartType;
#else 
  typedef ContainedGridPartType  GridPartType ;
#endif

typedef TestFunctionSpace FunctionSpaceType;
typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;

#if defined USE_BLOCKVECTORFUNCTION
typedef BlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#elif defined USE_VECTORFUNCTION
typedef ManagedDiscreteFunction
  < VectorDiscreteFunction
    < DiscreteFunctionSpaceType,
      Fem :: DynamicVector< FunctionSpaceType :: RangeFieldType > > >
  DiscreteFunctionType;
#elif defined USE_ATTACHEDFUNCTION
typedef AttachedDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#elif defined USE_BLOCKVECTORDISCRETEFUNCTION
typedef Dune::Fem::ReferenceBlockVector< double, DiscreteFunctionSpaceType::localBlockSize > 
  BlockVectorType;
typedef Dune::Fem::BlockVectorDiscreteFunction< DiscreteFunctionSpaceType, BlockVectorType > 
  DiscreteFunctionType;
#elif defined USE_COMBINEDFUNCTION
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
  ContainedDiscreteFunctionType;
typedef Dune::CombinedDiscreteFunction< ContainedDiscreteFunctionType, DIMRANGE >
  DiscreteFunctionType;
#else
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#endif

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

// dummy class for easier use of L2 Projection 
template< class Domain, class Range >
class DGL2Projection : 
  public L2Projection< typename Range :: DomainFieldType ,
                       typename Range :: RangeFieldType , 
                       Domain , Range > 
{
  typedef L2Projection< typename Range :: DomainFieldType ,
                       typename Range :: RangeFieldType , 
                       Domain , Range >  BaseType; 
public:
  DGL2Projection() : BaseType() {} 
};

// main program 
int main(int argc, char ** argv) 
{
  MPIManager :: initialize( argc, argv );
  try
  {
    MyGridType &grid = TestGrid :: grid();
    const int step = TestGrid :: refineStepsForHalf();
    grid.globalRefine( 2*step );

    GridPartType gridPart( grid );

    // add check for grid width 
    std::cout << "Grid width: " 
      << GridWidth :: calcGridWidth( gridPart ) << std::endl; 

    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    ExactSolutionType exactSolution;
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    solution.clear();
  
    // perform the L2Projection
    DGL2Projection< ExactSolutionType, DiscreteFunctionType > dgl2; 
    dgl2( exactSolution, solution );

    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
