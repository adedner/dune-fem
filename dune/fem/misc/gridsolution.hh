#ifndef DUNE_FEM_GRIDSOLUTION_HH
#define DUNE_FEM_GRIDSOLUTION_HH

#include <dune/common/exceptions.hh>
#include <dune/grid/utility/hierarchicsearch.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/io/file/datawriter.hh>

#if HAVE_DUNE_SPGRID 
#include <dune/grid/spgrid.hh>
#endif

namespace Dune {

namespace Fem {
  
// GridSolution
//-------------

/** \class GridSolution
 *  \brief creates a function with evaluate method from a check point
 *
 *  \tparam GridImp Grid type
 *  \tparam DiscreteFunctionImp Discrete function type
 */
template <class GridImp, class DiscreteFunctionImp > 
class GridSolution
{
  typedef GridImp  GridType;
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType ;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType ;
  typedef typename DiscreteFunctionSpaceType :: GridPartType   GridPartType ;
  typedef typename GridPartType :: IndexSetType IndexSetType;
  typedef CheckPointer< GridType >   CheckPointerType;

  typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType;
  typedef typename GridType :: template Codim<0> :: Entity        EntityType;

  typedef HierarchicSearch< GridType, IndexSetType > HierarchicSearchType;

  typedef tuple< DiscreteFunctionType* > IOTupleType;

  GridType* grid_;
  GridPtr< GridType > gridPtr_;
  GridPartType gridPart_;
  DiscreteFunctionSpaceType space_;
  DiscreteFunctionType discreteFunction_;
  IOTupleType data_;
  HierarchicSearchType hierarchicSearch_;

public:  
  GridType& grid() { assert( grid_ ); return *grid_ ; }
  const GridType& grid() const { assert( grid_ ); return *grid_ ; }

  //! Constructor
  explicit GridSolution(const std::string checkPointFile,
                        const int rank = -1 ) :
    grid_( CheckPointerType :: restoreGrid( checkPointFile, rank ) ), 
    gridPtr_(),
    gridPart_( grid() ),
    space_( gridPart_ ),
    discreteFunction_("grid-sol", space_),
    data_( &discreteFunction_ ),
    hierarchicSearch_( grid(), gridPart_.indexSet() )
  {
    // store grid pointer
    gridPtr_ = grid_;
    // restore data from checkpoint
    CheckPointerType :: restoreData( grid(), data_, checkPointFile, rank );
  }

  /** \brief evaluates in a given space-time point
   *  \param[in] x Point in global coordinates
   *  \param[in] time Time
   *  \param[out] result The value of the discrete function in space-time point \f[ (x,time)\f]
   *
   *  \tparam PointType The point type
   */
  void evaluate(const DomainType& x, const double time, RangeType& result) const 
  {
    evaluate(x, result );
  }

  /** \brief evaluates in a given space point
   *  \param[in] x Point in global coordinates
   *  \param[out] result The value of the discrete function in space point @a x
   *
   *  \tparam PointType The point type
   */
  void evaluate(const DomainType& x, RangeType& result) const 
  {
    // search entity 
    EntityPointerType ep = hierarchicSearch_.findEntity( x );
    const EntityType& entity = *ep ;

    const DomainType local = entity.geometry().local( x );

    // evaluate discrete function 
    discreteFunction_.localFunction( entity ).evaluate( local, result );
  }

  //! \brief writes a discrete function 
  static void writeDiscreteFunction(const GridType& grid, 
                                    const DiscreteFunctionType& discreteFunction,
                                    const double time ) 
  {
    typedef tuple< const DiscreteFunctionType* > OutputTuple;
    OutputTuple data( &discreteFunction );
    CheckPointerType :: writeSingleCheckPoint( grid, data, time, false );
  }
};

///////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////
/** \class GridSolution
 *  \brief creates a function with evaluate method from a check point
 *
 *  \tparam GridImp Grid type
 *  \tparam DiscreteFunctionImp Discrete function type
 */
template <class GridImp, class DiscreteFunctionImp > 
class GridSolutionVector
{
  typedef GridImp  GridType;
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType :: DomainType  DomainType;
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType :: RangeType   RangeType;

  typedef GridSolution< GridType, DiscreteFunctionType > GridSolutionType ;

  template <class DomainType, class Grid>
  struct CheckDomain
  {
    static bool isInside(const DomainType& x, const Grid& grid ) 
    {
      typedef typename Grid :: LevelGridView MacroView ;
      typedef typename MacroView :: template Codim< 0 > :: Iterator Iterator ;
      typedef typename Iterator :: Entity Entity;
      const MacroView& macroView = grid.levelView( 0 );
      const Iterator end = macroView.template end< 0 > ();
      for( Iterator it = macroView.template begin< 0 > (); it != end; ++it )
      {
        const Entity& entity = * it ;
        // check that corners are within the reference element of the given type
        const GenericReferenceElement< typename Grid::ctype, Grid::dimensionworld > &refElement
             = GenericReferenceElements< typename Grid::ctype, Grid::dimensionworld >::general( entity.type() );
        
        typedef typename Entity :: Geometry Geometry;
        const Geometry& geo = entity.geometry();

        if( refElement.checkInside( geo.local( geo.center() ) ) )
          return true ;
      }
      return false ;
    }
  };

#if HAVE_DUNE_SPGRID 
  template <class DomainType, class ct, int dim, SPRefinementStrategy strategy , class Comm>
  struct CheckDomain< DomainType, SPGrid< ct, dim, strategy, Comm > > 
  {
    typedef SPGrid< ct, dim, strategy, Comm >  Grid;
    static bool isInside(const DomainType& x, const Grid& grid ) 
    {
      return grid.domain().contains( x );
    }
  };
#endif


  const int numProcs_;
  std::vector< GridSolutionType* > solutions_; 

  const int numProcs(const std::string& checkPointFile) const 
  {
    int numProc = MPIManager :: size();
    readParameter(checkPointFile,"NumberProcessors",numProc,Parameter::verbose ());
    return numProc;
  }
public:  
  //! Constructor
  explicit GridSolutionVector(const std::string checkPointFile) :
    numProcs_( numProcs( checkPointFile ) ),
    solutions_( numProcs_, (GridSolutionType *) 0 )
  {
    for(int p=0; p<numProcs_; ++p)
    {
      std::cout << "Reading Grid " << p << " from checkpoint" << std::endl;
      solutions_[ p ] = new GridSolutionType( checkPointFile, p );
    }
  }

  ~GridSolutionVector() 
  {
    for(int p=0; p<numProcs_; ++p)
    {
      delete solutions_[ p ];
      solutions_[ p ] = 0;
    }
  }

  /** \brief evaluates in a given space-time point
   *  \param[in] x Point in global coordinates
   *  \param[in] time Time
   *  \param[out] result The value of the discrete function in space-time point \f[ (x,time)\f]
   *
   *  \tparam PointType The point type
   */
  void evaluate(const DomainType& x, const double time, RangeType& result) const 
  {
    evaluate(x, result );
  }

  /** \brief evaluates in a given space point
   *  \param[in] x Point in global coordinates
   *  \param[out] result The value of the discrete function in space point @a x
   *
   *  \tparam PointType The point type
   */
  void evaluate(const DomainType& x, RangeType& result) const 
  {
    for(int p=0; p<numProcs_; ++p)
    {
      assert( solutions_[ p ] );
      const GridSolutionType& gridSolution = *(solutions_[ p ]);
      if( isInDomain( x, gridSolution.grid() ) )
      {
        gridSolution.evaluate( x, result );
        return ;
      }
    }

    std::cerr << "Point in Grid not found" << std::endl;
    assert( false );
    abort();
  }

  bool isInDomain( const DomainType& x, const GridType& grid ) const 
  {
    return CheckDomain<DomainType, GridType> :: isInside( x, grid );
  }

  //! \brief writes a discrete function 
  static void writeDiscreteFunction(const GridType& grid, 
                                    const DiscreteFunctionType& discreteFunction,
                                    const double time ) 
  {
    GridSolutionType :: writeDiscreteFunction( grid, discreteFunction, time );
  }
};

} // end namespace Fem 
} // end namespace Dune 
#endif
