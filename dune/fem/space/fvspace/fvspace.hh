#ifndef DUNE_FVSPACE_HH
#define DUNE_FVSPACE_HH

#include <map>

#include <dune/grid/common/grid.hh>

#include <dune/fem/space/common/discretefunctionspace.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/basefunctions/basefunctionstorage.hh>
#include <dune/fem/space/basefunctions/basefunctionsets.hh>
#include <dune/fem/space/basefunctions/basefunctionproxy.hh>

#include <dune/fem/space/fvspace/fvspacebasefunctions.hh>

#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include <dune/fem/space/common/defaultcommhandler.hh>

// * Note: the dofmanager could be removed from the space altogether now.
// (Maybe this wouldn't be a clever move, though. In my view of a perfect Dune,
// there would be one DofManager per space and the DiscreteFunctions wouldn't
// need to fiddle with the DofMapper anymore...

namespace Dune
{

  // Forward declarations
  template <class FunctionSpaceImp, class GridPartImp, int polOrd,
            template<class> class BaseFunctionStorageImp = CachingStorage >
  class FiniteVolumeSpace;

  template <class FunctionSpaceImp,class GridPartImp, int polOrd,
            template <class> class BaseFunctionStorageImp > 
  struct FiniteVolumeSpaceTraits 
  {   
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef GridPartImp GridPartType;

    typedef typename GridPartType::GridType GridType;
    typedef typename GridPartType::IndexSetType IndexSetType;
    typedef typename GridPartType::template Codim<0>::IteratorType IteratorType;
    enum { dimRange = FunctionSpaceType :: dimRange };

    // dimension of local coordinates 
    enum { dimLocal = GridType :: dimension };

    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef FiniteVolumeSpace<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceType;
    
    // convert function space to local function space 
    typedef typename ToLocalFunctionSpace< FunctionSpaceImp, dimLocal > :: Type
      BaseFunctionSpaceType ; 

    typedef VectorialBaseFunctionSet<BaseFunctionSpaceType, BaseFunctionStorageImp > BaseFunctionSetImp;

    typedef SimpleBaseFunctionProxy<BaseFunctionSetImp> BaseFunctionSetType;

    enum { localBlockSize = dimRange };

    // block mapper 
    typedef CodimensionMapper< GridPartType, 0 > BlockMapperType;

    // type of mapper for block vector functions 
    typedef NonBlockMapper< BlockMapperType, localBlockSize > MapperType;

    
    /** \brief defines type of data handle for communication 
        for this type of space.
    */
    template< class DiscreteFunction,
              class Operation = DFCommunicationOperation :: Copy >
    struct CommDataHandle
    {
      //! type of data handle 
      typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      //! type of operation to perform on scatter 
      typedef Operation OperationType;
    };
  };
  //
  //  --FiniteVolumeSpace
  //

  /** @addtogroup FVDFSpace

   Provides access to base function set for different element 
   type in one grid and size of functionspace 
   and map from local to global dof number

   \note This space can only be used with a special set of index sets.
   If you want to use the FiniteVolumeSpace with an index set only
   supportting the index set interface, then use the IndexSetWrapper
   class which will add the needed functionalty.

   \note For adaptive calculations one have to use Index Sets that are
   capable for adaptation, i.e. the method adaptive returns true, see 
   AdaptiveLeafIndexSet. 
   @{
  **/

  /** @brief 
      Finite Volume Function Space 
      **/
  template<class FunctionSpaceImp, class GridPartImp, int polOrd, 
           template <class> class BaseFunctionStorageImp >
  class FiniteVolumeSpace : 
    public DiscreteFunctionSpaceDefault
  <
    FiniteVolumeSpaceTraits<FunctionSpaceImp, GridPartImp, 
                            polOrd, BaseFunctionStorageImp > 
  >
  {
 public:
    typedef typename GridPartImp::GridType GridType;

    typedef FiniteVolumeSpace< 
          FunctionSpaceImp, GridPartImp, polOrd , BaseFunctionStorageImp
      > FiniteVolumeSpaceType;

    //! type of this pointer 
    typedef FiniteVolumeSpaceType ThisType; 
 
    //! my Traits 
    typedef FiniteVolumeSpaceTraits<
      FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp
      > Traits;

    typedef DiscreteFunctionSpaceDefault<Traits> DefaultType;
  
    /** type of base function set implementation  */
    typedef typename Traits::BaseFunctionSetImp BaseFunctionSetImp;

    typedef typename Traits::BaseFunctionSetType BaseFunctionSetType;

    typedef typename Traits::IndexSetType IndexSetType;

    typedef typename Traits::GridPartType GridPartType;
    
    typedef typename Traits::IteratorType IteratorType;

    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    //! type of id 
    typedef int IdentifierType;
    
    //! id is neighbor of the beast
    static const IdentifierType id = 665;

    typedef typename Traits :: MapperType MapperType; 

    //! block mapper 
    typedef typename Traits :: BlockMapperType BlockMapperType; 

 protected:  
    //! mapper factory 
    typedef CodimensionMapperSingletonFactory< GridPartType, 0 > BlockMapperSingletonFactoryType;

    //! singleton list of mappers 
    typedef SingletonList
      < typename BlockMapperSingletonFactoryType::Key, BlockMapperType, BlockMapperSingletonFactoryType >
      BlockMapperProviderType;

    //! scalar basefunction space type 
    typedef typename Traits::BaseFunctionSpaceType BaseFunctionSpaceType;

  public:
    //! type of base function factory 
    typedef FVBaseFunctionFactory<typename BaseFunctionSpaceType ::
      ScalarFunctionSpaceType, polOrd> ScalarFactoryType;

    //! type of singleton factory 
    typedef BaseFunctionSetSingletonFactory<GeometryType,BaseFunctionSetImp,
                ScalarFactoryType> SingletonFactoryType;

    //! type of singleton list  
    typedef SingletonList< GeometryType, BaseFunctionSetImp,
            SingletonFactoryType > SingletonProviderType;

    //! default communication interface 
    static const InterfaceType defaultInterface = InteriorBorder_All_Interface;

    //! default communication direction 
    static const CommunicationDirection defaultDirection =  ForwardCommunication;

    //! remember polynomial order 
    enum { polynomialOrder =  polOrd };

    //! Constructor generating Finite Volume Space 
    inline explicit FiniteVolumeSpace(GridPartType & g,
          const InterfaceType commInterface = defaultInterface ,
          const CommunicationDirection commDirection = defaultDirection );

    //! Desctructor 
    ~FiniteVolumeSpace (); 

    /** \copydoc Dune::DiscreteFunctionSpaceInterface::contains */
    inline bool contains(const int codim) const
    {   
      return blockMapper().contains( codim );
    }
    
    /** \copydoc Dune::DiscreteFunctionSpaceInterface::continuous */
    bool continuous() const { return (polynomialOrder == 0) ? false : true; }
 
    //! return type of this function space 
    DFSpaceIdentifier type () const;

    //! returns polynomial order
    int order() const { return polynomialOrder; }

    //! provide the access to the base function set for a given entity
    template <class EntityType>
    const BaseFunctionSetType 
    baseFunctionSet ( const EntityType &en ) const;

    //! Get base function set for a given id of geom type (mainly used by
    //! CombinedSpace) 
    const BaseFunctionSetType
    baseFunctionSet (const GeometryType& geoType) const;

    //! get dimension of value 
    int dimensionOfValue () const;

    //! Return the dof mapper of the space
    MapperType& mapper() const;

    //! Return the dof mapper of the space
    BlockMapperType& blockMapper() const;

  protected:
    //! create functions space
    void makeFunctionSpace (GridPartType& gridPart); 
  
  protected:
    //! type of corresponding map of base function sets
    typedef std::map< const GeometryType, const BaseFunctionSetImp *> BaseFunctionMapType;

    //! the corresponding map of base function sets
    mutable BaseFunctionMapType baseFuncSet_;

  private:
    //! mapper for block vector functions 
    BlockMapperType& blockMapper_;
    //! the corresponding FiniteVolumeMapper 
    mutable MapperType mapper_; 
  }; // end class FiniteVolumeSpace


/** \brief This is a restriction/prolongation operator for FV data of order zero.
 */
template <class DiscreteFunctionImp,
          class FunctionSpaceImp, 
          class GridPartImp, 
          template <class> class StorageImp> 
class RestrictProlongDefaultImplementation< 
  DiscreteFunctionImp,
  FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, 0,StorageImp> 
  > :
public RestrictProlongPieceWiseConstantData<DiscreteFunctionImp> 
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef RestrictProlongPieceWiseConstantData<DiscreteFunctionImp> BaseType; 
public:  
  //! Constructor
  RestrictProlongDefaultImplementation( DiscreteFunctionType & df ) : 
    BaseType ( df ) 
  {
  }
};

/** @} **/  
} // end namespace Dune

// contains the implementation of FiniteVolumeSpace
#include "fvspace_inline.hh"
#endif
