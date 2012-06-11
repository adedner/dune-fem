#ifndef DUNE_ADAPTMANAGER_HH
#define DUNE_ADAPTMANAGER_HH

//- system includes 
#include <string>

//- local includes 
#include <dune/common/timer.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/operator/common/objpointer.hh>

#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/loadbalancer.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/storage/singletonlist.hh>

#include <dune/fem/space/common/adaptcallbackhandle.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/threadmanager.hh>


namespace Dune
{

namespace Fem 
{

/** @addtogroup Adaptation Adaptation 
    Here the interfaces and algorithms for adapatation of a grid are
    described and implemented. 

    The strategy for restrict/prolong of data is chosen through the
    parameter
     \parametername \c fem.adaptation.method \n
                    method used for 
                    prolongation/restriction of data during grid
                    refinement;
                    defaults to generic \n
    Values are:
     - 0 == none: no adaptation is performed
     - 1 == generic: works on all grids 
     - 2 == callback: only AlbertaGrid and ALUGrid
     .
     
    \remarks 
    The Interface of Adaptation Manager is described by the class
    AdaptationManagerInterface.
    @{
 **/

/** \brief AdaptationManagerInterface class. 
 
 This Class is the result of a combination of different
 AdaptationOperators. It is the same principle as with Mapping. 
*/ 
class AdaptationManagerInterface : virtual public LoadBalancerInterface 
{
public:
  //! \brief default constructor 
  AdaptationManagerInterface () : am_ (0) {}

  //! \brief destructor 
  virtual ~AdaptationManagerInterface () {}
  
  /** \brief on call of this method the internal adaptation operator is
      called. 
  */
  virtual void adapt ()  
  {
    //std::cout << "called AdaptationManagerInterface::adapt()" << std::endl;
    if(am_) am_->adapt();  
    else 
    {
      std::cerr << "WARNING: AdaptationManagerInterface::adapt: no adaptation manager assigned! \n";
    }
  }

  /** \brief returns true if adaptation manager as adaptation method different to NONE
     \return \b true if adaptation method is not NONE, \b false otherwise
  */
  virtual bool adaptive () const  
  { 
    return (am_) ? (am_->adaptive()) : false; 
  } 

  /** \brief returns name of adaptation method 
     \return name of adaptation method 
  */
  virtual const char * methodName() const 
  {
    return (am_) ? (am_->methodName()) : "unknown method";
  }
    
  /** \brief Assignment operator, pointer to adaptation manager is stored 
      \return reference to this (i.e. *this) 
  */
  AdaptationManagerInterface & operator = (const AdaptationManagerInterface & am)
  {
      /** \todo This const-casting seems strange to me! */
    am_ = const_cast<AdaptationManagerInterface *> (&am);
    return (*this);
  }

  /** @copydoc LoadBalancerInterface::loadBalance */
  virtual bool loadBalance () 
  { 
    return (am_) ? (am_->loadBalance()) : false; 
  }

  /** @copydoc LoadBalancerInterface::balanceCounter */
  virtual int balanceCounter () const 
  { 
    return (am_) ? (am_->balanceCounter()) : 0; 
  }

  /** \brief time that last adaptation cycle took */
  virtual double adaptationTime () const 
  {
    return 0.0;
  }
 
private: 
  //! pointer to adaptation manager 
  AdaptationManagerInterface* am_;
};

/** \brief AdaptationMethod is a simple adaptation method reader class. */
template <class GridType>
class AdaptationMethod : virtual public AdaptationManagerInterface 
{  
public:  
  //! type of adaptation method 
  enum AdaptationMethodType { none = 0, //!< no adaptation is performed 
                              generic = 1, //!< a generic restriction and prolongation algorithm is used  
                              callback = 2 //!< the callback mechanism from AlbertaGrid and ALUGrid is used 
                            };
  
public:
  /** \brief constructor of AdaptationMethod 
     The following optional parameters are used 
        # 0 == none, 1 == generic, 2 == call back (only AlbertaGrid and ALUGrid)  
        AdaptationMethod: 1 # default value 
     \param grid Grid that adaptation method is read for  

  */   
  AdaptationMethod(const GridType & grid) 
    : adaptationMethod_(generic) 
  {
    const bool output = ( Parameter :: verbose() );
    int am = 1;
    const std::string methodNames [] = { "none", "generic", "callback" };
    am = Parameter::getEnum("fem.adaptation.method", methodNames, am);
    init(am,output);
  }
private:
  void init(int am,const bool output)
  {

    // chose adaptation method 
    if(am == 2) adaptationMethod_ = callback;
    else if(am == 1) adaptationMethod_ = generic;
    else adaptationMethod_ = none;

    // for structred grid adaptation is disabled 
    if( ! Capabilities::isLocallyAdaptive<GridType>::v ) 
    {
      adaptationMethod_ = none;
      if( output )
      {
        std::cerr << "WARNING: AdaptationMethod: adaptation disabled for structured grid! \n";
      }
    }
      
    if( output )
    {
      std::cout << "Created AdaptationMethod: adaptation method = " << methodName() << std::endl;
    }
  }
public:
  //! virtual destructor 
  virtual ~AdaptationMethod () {}

  /** \copydoc Dune::AdaptationManagerInterface::methodName */
  virtual const char * methodName() const 
  {
    switch (adaptationMethod_) {
      case generic: return "generic";
      case callback: return "callback";              
      case none: return "no adaptation";
      default:  return "unknown method";
    }
  }

  /** \copydoc Dune::AdaptationManagerInterface::adaptive */
  virtual bool adaptive () const { return adaptationMethod_ != none; }
  
protected:  
  //! method identifier 
  AdaptationMethodType adaptationMethod_;
};

/*! This could be seen as a hack, but is not
  With this method we define the class CombineRestProl which is only 
  for combining of local grid operations without using virtual methods.
  The interface of the two defined methods of the class (restrictLocal and
  prolongLocal) is given by the implementation (see below ) and 
  has to be the same for all local operators you want to combine.
*/
#define PARAM_CLASSNAME CombinedRestProl 
#define PARAM_FUNC_1 restrictLocal 
#define PARAM_FUNC_2 prolongLocal 
#define PARAM_FUNC_3 calcFatherChildWeight 
#define PARAM_FUNC_4 addToCommunicator 
#include <dune/fem/operator/common/combine.inc>

/*! \brief This class manages the adaptation process. 
 If the method adapt is called, then the grid is adapted and also 
 all the data belonging to the given dof manager will be rearranged 
 for data set where it is necessary to keep the data.
 */
template <class GridType, class RestProlOperatorImp >
class AdaptationManagerBase
: public AdaptationMethod< GridType >,
  public ObjPointerStorage 
{  
  typedef AdaptationMethod< GridType > BaseType;
  typedef typename BaseType :: AdaptationMethodType AdaptationMethodType; 
  
  template <class AdaptManager, class GridImp, bool isGoodGrid> 
  struct CallAdaptationMethod
  {
    template <class DofManagerImp, class RPOpImp>
    static void adapt(const AdaptManager& am, GridImp & grid, 
                      DofManagerImp& dm , RPOpImp& rpop,
                      AdaptationMethodType adaptMethod) 
    {
      // use generic adapt method 
      if( adaptMethod == BaseType :: generic ) 
      {
        am.template genericAdapt<All_Partition> ();
        return ;
      }
      
      // use grid call back adapt method 
      if( adaptMethod == BaseType :: callback ) 
      {
        // combine dof manager and restrict prolong operator 
        typedef RestrictProlongWrapper< GridImp, DofManagerType, RPOpImp > RPType;

        // create new handle 
        RPType restrictProlongHandle ( dm , rpop );

        // call grid adaptation 
        grid.adapt( restrictProlongHandle ); 
        return ;
      }
    }
  };
  
  template <class AdaptManager, class GridImp> 
  struct CallAdaptationMethod<AdaptManager,GridImp,false>
  {
    template <class DofManagerImp, class RPOpImp>
    static void adapt(const AdaptManager& am, GridImp & grid, 
                      DofManagerImp& dm , RPOpImp& rpop,
                      AdaptationMethodType adaptMethod) 
    {
      // use generic adapt method 
      if(adaptMethod != BaseType :: none ) 
      {
        // use partition type All_Partition, 
        // since we also need to iterate on ghosts 
        // for wasChanged information gathering
        am.template genericAdapt<All_Partition> ();
        return ;
      }
    }
  };
  
  //! type of this class 
  typedef AdaptationManagerBase<GridType,RestProlOperatorImp> ThisType;

  //! type of dof manager 
  typedef DofManager< GridType > DofManagerType; 

public:
  typedef typename GridType :: Traits :: LocalIdSet LocalIdSet;

  /** \brief constructor of AdaptationManagerBase 
     The following optional parameter can be used
        # 0 == none, 1 == generic, 2 == call back (only AlbertaGrid and ALUGrid)  
        AdaptationMethod: 1 # default value 
     \param grid Grid that adaptation is done for 
     \param rpOp restriction and prlongation operator that describes how the 
      user data is projected to other grid levels
        the following two lines:
  */   
  AdaptationManagerBase ( GridType &grid, RestProlOperatorImp &rpOp )
  : BaseType( grid ),
    grid_( grid ),
    dm_( DofManagerType::instance( grid_ ) ),
    rpOp_( rpOp ),
    adaptTime_( 0.0 ),
    wasChanged_( false )
  {}

  //! destructor 
  virtual ~AdaptationManagerBase () {}

  //! no public method, but has to be public, because all AdaptationManagers 
  //! must be able to call this method and the template parameters are
  //! allways different 
  RestProlOperatorImp & getRestProlOp ()  
  {
    return rpOp_;
  }
  
  /** \brief 
     according to adaption method parameter 
     the adaption procedure is done, 
     0 == no adaptation
     1 == generic adaption 
     2 == grid call back adaptation (only in AlbertaGrid and ALUGrid)
  */
  virtual void adapt ()
  {
    // make sure this is only called in single thread mode 
    assert( Fem :: ThreadManager :: singleThreadMode() );

    // get stopwatch 
    Timer timer; 

    const bool supportsCallback = Capabilities :: supportsCallbackAdaptation< GridType > :: v;
    CallAdaptationMethod< ThisType, GridType, supportsCallback >
      :: adapt(*this,grid_,dm_,rpOp_,this->adaptationMethod_);
    
    // take time 
    adaptTime_ = timer.elapsed();
  }

  //! \brief default load balancing method which does nothing
  virtual bool loadBalance () 
  { 
    return false; 
  }

  //! default load balancing counter is zero 
  virtual int balanceCounter () const 
  { 
    return 0; 
  }

  /** \copydoc Dune::AdaptationManagerInterface::adaptationTime */
  virtual double adaptationTime() const 
  {
    return adaptTime_;
  }

protected:
  static DofManagerType& getDofManager(const GridType& grid) 
  {
    return DofManagerType :: instance( grid );
  }

private:  
  /** \brief generic adaptation procedure
     adapt defines the grid walkthrough before and after grid adaptation.
     Note that the LocalOperator can be an combined Operator 
     Domain and Range are defined through class Operator
  */
  template <PartitionIteratorType pitype>
  void genericAdapt () const 
  {
    // call pre-adapt, returns true if at least 
    // one element is marked for coarsening 
    bool restr = grid_.preAdapt();  

    // get macro grid view  
    typedef typename GridType::template Partition< All_Partition > :: LevelGridView  MacroGridView ;
    typedef typename MacroGridView :: template Codim<0>::
            template Partition<pitype> :: Iterator MacroIterator;

    // reset flag 
    wasChanged_ = false ;

    if(restr)
    {
      // get macro grid view 
      MacroGridView macroView = grid_.levelView( 0 );

      // make a hierarchical to insert all elements 
      // that are father of elements that might be coarsened 

      {
        // get macro iterator 
        MacroIterator endit  = macroView.template end<0,pitype>  ();
        for(MacroIterator it = macroView.template begin<0,pitype>();
            it != endit; ++it )
        {
          hierarchicRestrict( *it , dm_.indexSetRestrictProlongNoResize() );
        }
      }

      // if at least one element was found for restriction 
      if( wasChanged_ )
      {
        // now resize memory 
        dm_.resizeForRestrict();

        // now project all data to fathers 
        {
          // get macro iterator 
          MacroIterator endit  = macroView.template end<0,pitype>  ();
          for(MacroIterator it = macroView.template begin<0,pitype>();
              it != endit; ++it )
          {
            hierarchicRestrict( *it , rpOp_ );
          }
        }
      }
    }
    
    // adapt grid due to preset markers
    // returns true if at least one element was refined 
    const bool refined = grid_.adapt();

    if( refined )
    {
      // get macro grid view 
      MacroGridView macroView = grid_.levelView( 0 );

      // resizes the index sets (insert all new indices) 
      // and resizes the memory
      dm_.resize();

      // make run through grid to project data 
      MacroIterator endit  = macroView.template end<0,pitype>  ();
      for(MacroIterator it = macroView.template begin<0,pitype>();
          it != endit; ++it )
      {
        hierarchicProlong( *it , rpOp_ );
      }
    }

    // notifyGlobalChange make wasChanged equal on all cores 
    if( dm_.notifyGlobalChange( wasChanged_ ) )
    {
      // compress index sets and data 
      // this will increase the sequence counter
      dm_.compress();
    }

    // do cleanup 
    grid_.postAdapt();
  }
  
private:
  //! make hierarchic walk trough for restriction 
  template <class EntityType, class RestrictOperatorType  >
  bool hierarchicRestrict ( const EntityType& entity, RestrictOperatorType & restop ) const 
  {
    if( ! entity.isLeaf() )
    {
      // true means we are going to restrict data 
      bool doRestrict = true;

      // check partition type 
      const bool isGhost = entity.partitionType() == GhostEntity ;
      
      // if the children have children then we have to go deeper 
      const int childLevel = entity.level() + 1;
      typedef typename EntityType::HierarchicIterator HierarchicIterator; 
      
      // check all children first 
      {
        const HierarchicIterator endit = entity.hend( childLevel );
        for(HierarchicIterator it = entity.hbegin( childLevel ); it != endit; ++it)
        {
          doRestrict &= hierarchicRestrict( *it , restop );
        }
      }

      // if doRestrict is still true, restrict data 
      if(doRestrict)
      {
        // we did at least one restriction 
        wasChanged_ = true; 

        // do not restrict the solution on ghosts, this will 
        // fail, but we still need the wasChanged info, so simply 
        // calling hierarchicRestrict on interior won't work either
        if( ! isGhost ) 
        {
          // true for first child, otherwise false 
          bool initialize = true;
          const HierarchicIterator endit = entity.hend( childLevel );
          for(HierarchicIterator it = entity.hbegin( childLevel ); it != endit; ++it)
          {
            // restrict solution 
            restop.restrictLocal( entity, *it , initialize);
            // reset initialize flag 
            initialize = false;
          }
        }
      }
    }

    // if all children return mightBeCoarsened,
    // then doRestrict on father remains true 
    return entity.mightVanish();
  }

  template <class EntityType, class ProlongOperatorType >
  void hierarchicProlong ( const EntityType &entity, ProlongOperatorType & prolop ) const 
  {
    typedef typename EntityType::HierarchicIterator HierarchicIterator; 
    typedef typename GridType :: template Codim<0> :: EntityPointer EntityPointerType; 
    
    // NOTE: initialize not working here
    // because we call hierarchically  
    
    // first call on this element 
    bool initialize = true;

    // check partition type 
    const bool isGhost = entity.partitionType() == GhostEntity ;
    
    const int maxLevel = grid_.maxLevel();
    const HierarchicIterator endit = entity.hend( maxLevel );
    for( HierarchicIterator it = entity.hbegin( maxLevel ); it != endit; ++it )
    {
      // should only get here on non-leaf entities 
      assert( !entity.isLeaf() );

      const EntityType & son = *it; 
      if( son.isNew() )
      {
        // the grid was obviously changed if we get here
        wasChanged_ = true ;

        // do not prolong the solution on ghosts, this will 
        // fail, but we still need the wasChanged info, so simply 
        // calling hierarchicRestrict on interior won't work either
        if( ! isGhost ) 
        {
          // get father entity 
          EntityPointerType vati = son.father();

          // do prolongation 
          prolop.prolongLocal( *vati , son , initialize ); 
          initialize = false;
        }
      }
    }
  }
 
protected:  
  //! corresponding grid 
  GridType &grid_;

  //! DofManager corresponding to grid
  DofManagerType &dm_;
  
  //! Restriction and Prolongation Operator 
  RestProlOperatorImp &rpOp_;

  //! time that adaptation took 
  double adaptTime_;

  //! flag for restriction 
  mutable bool wasChanged_; 
};

//! factory class to create adaptation manager reference counter 
template <class KeyType, class ObjectType>
struct AdaptationManagerReferenceFactory
{
  static ObjectType* createObject(const KeyType& key)
  {
    return new ObjectType(0);
  }
  static void deleteObject(ObjectType* obj)
  {
    delete obj;
  }
};

/*! \brief This class manages the adaptation process including a load
  balancing after the adaptation step. This class is created by the
  AdaptationManager for each grid instance. See AdaptationManager for
  details. 
*/
template <class GridType, class RestProlOperatorImp>
class AdaptationManager :
  public AdaptationManagerBase<GridType,RestProlOperatorImp> ,
  public LoadBalancer<GridType> 
{  
  // type of key 
  typedef const GridType* KeyType;
  // object type 
  typedef size_t ObjectType;
  // type of factory 
  typedef AdaptationManagerReferenceFactory<KeyType, ObjectType>  FactoryType;

  // type of singleton list 
  typedef SingletonList< KeyType, ObjectType, FactoryType > ProviderType;

  typedef AdaptationManagerBase<GridType,RestProlOperatorImp> BaseType; 
  typedef LoadBalancer<GridType> Base2Type;

  mutable CommunicationManagerList commList_;
  double balanceTime_ ;

  // reference counter to ensure only one instance per grid exists 
  ObjectType& referenceCounter_;

  // do not copy 
  AdaptationManager(const AdaptationManager&);

public:  
  /** \brief constructor of AdaptationManager 
     The following optional parameters from the Dune::Parameter class are used
        # 0 == none, 1 == generic, 2 == call back (only AlbertaGrid and ALUGrid)  
        fem.adaptation.method: 1 # default value 
     \param grid Grid that adaptation is done for 
     \param rpOp restriction and prlongation operator that describes how the 
      user data is projected to other grid levels
     \param balanceCounter start counter for balance cycle (default = 0)   
  */   
  AdaptationManager(GridType & grid, 
                    RestProlOperatorImp & rpOp, 
                    int balanceCounter = 0 )
    : BaseType(grid,rpOp) 
    , Base2Type( grid , rpOp , balanceCounter )
    , commList_(rpOp)
    , balanceTime_( 0.0 )
    , referenceCounter_( ProviderType :: getObject( &grid ) )
  {
    ++ referenceCounter_;
    if( referenceCounter_ > 1 )
    {
      DUNE_THROW(InvalidStateException,"Only one instance of AdaptationManager allowed per grid instance");
    }
  }

  //! destructor decreasing reference counter 
  ~AdaptationManager() 
  {
    -- referenceCounter_;
    ProviderType :: removeObject( referenceCounter_ );
  }

  /** @copydoc LoadBalancerInterface::loadBalance */
  virtual bool loadBalance () 
  {
    // call load balance 
    return Base2Type :: loadBalance();
  }

  /** @copydoc LoadBalancerInterface::balanceCounter */ 
  virtual int balanceCounter () const 
  { 
    return Base2Type :: balanceCounter ();
  }

  /** @copydoc LoadBalancerInterface::loadBalanceTime */
  virtual double loadBalanceTime() const 
  {
    return balanceTime_;
  }
  
  /** @copydoc AdaptationManagerInterface::adapt */ 
  virtual void adapt () 
  {
    // adapt grid 
    BaseType :: adapt ();

    // if adaptation is enabled 
    if( this->adaptive() )
    {
      // get stopwatch 
      Timer timer; 
    
      // do load balancing 
      loadBalance ();

      // exchange all modified data 
      // this also rebuilds the dependecy cache of the 
      // cached communication manager if used 
      commList_.exchange();

      // get time  
      this->balanceTime_ = timer.elapsed();
    }
  }
};

/** \brief A class with one static method apply to globally refine a grid.
    All index sets are adapted to the new grid and the 
    managed dof storage is expanded - but no prolongation or
    restriction of data is performed.
*/
struct GlobalRefine 
{
  /** \brief apply global refinement and also adjust index sets and 
      managed dof storage. However, user data stored before is lost. 
      \param grid Grid that is globally refined 
      \param step refinement steps that are applied 
  */
  template <class GridType>
  static void apply(GridType& grid, const int step) 
  {
    typedef DofManager< GridType > DofManagerType; 
    DofManagerType& dm = DofManagerType :: instance(grid);
    grid.globalRefine(step);
    grid.loadBalance();
    dm.resize();
    dm.compress();
  }
};
/** \brief A class with one static method apply for invoking the local
    adaptation procedure on a given grid instance.
    All index sets are adapted to the new grid and the 
    managed dof storage is expanded - but no prolongation or
    restriction of data is performed.
*/
struct LocalRefine 
{
  /** \brief apply local refinement and also adjust index sets and 
      managed dof storage. However, user data stored before is lost. 
      \param grid Grid that is globally refined 
  */
  template <class GridType>
  static void apply(GridType& grid) 
  {
    typedef DofManager< GridType > DofManagerType; 
    DofManagerType& dm = DofManagerType :: instance(grid);
    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();
    grid.loadBalance();
    dm.resize();
    dm.compress();
  }
};

/** @} end documentation group */

} // namespace Fem

// #if DUNE_FEM_COMPATIBILITY  
// put this in next version 1.4 

using Fem :: AdaptationMethod ;
using Fem :: AdaptationManager ;
using Fem :: LocalRefine ;
using Fem :: GlobalRefine ;

// #endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune 

#endif // #ifndef DUNE_ADAPTMANAGER_HH
