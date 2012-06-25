#ifndef DUNE_LOCALFUNCTIONWRAPPER_HH
#define DUNE_LOCALFUNCTIONWRAPPER_HH

//-s system includes 
#include <cassert>

//- Dune includes 
#include <dune/fem/storage/objectstack.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{
  namespace Fem 
  {

  template< class DFTraits >
  class DiscreteFunctionDefault;

  template< class LocalFunctionStorage >
  class LocalFunctionWrapper;



  template< class LocalFunctionStorage >
  struct LocalFunctionWrapperTraits
  {
    typedef LocalFunctionStorage LocalFunctionStorageType;
    typedef typename LocalFunctionStorageType :: ObjectType
      LocalFunctionImpType;
    typedef typename LocalFunctionImpType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    typedef LocalFunctionWrapper< LocalFunctionStorageType >
      LocalFunctionUserType;
  };



  //**************************************************************************
  //
  //  --LocalFunctionWrapper 
  //
  //**************************************************************************
  //! Manages the getting and deleting of local function pointers and 
  //! acts like a local functions 
  template < class LocalFunctionStorage >
  class LocalFunctionWrapper
  : public LocalFunction< LocalFunctionWrapperTraits< LocalFunctionStorage > >
  {
  public:
    //! type of the traits
    typedef LocalFunctionWrapperTraits< LocalFunctionStorage > Traits;

    //! type of the local function storage (usually LocalFunctionStack, see below)
    typedef typename Traits :: LocalFunctionStorageType
      LocalFunctionStorageType;

    //! type of wrapped local function implementation 
    typedef typename Traits :: LocalFunctionImpType LocalFunctionImpType;

    //! type of discrete function space 
    typedef typename Traits :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

  private:
    typedef LocalFunctionWrapper< LocalFunctionStorageType > ThisType;
    typedef LocalFunction< Traits > BaseType;

    friend class EngineWrapper< LocalFunctionImpType, ThisType >;
    
  private:
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;

  public:
    //! type of the entity, the local function lives on is given by the space
    typedef typename DiscreteFunctionSpaceType :: EntityType  EntityType;

    enum {dimRange  = DiscreteFunctionSpaceType :: dimRange };
    enum {dimDomain = DiscreteFunctionSpaceType :: dimDomain };

  protected:
    // type of stack entry 
    typedef typename LocalFunctionStorageType :: PointerType
      LocalFunctionImpPtrType;
   
  protected:
    // pair storing pointer to local function and poiner to ref-counter 
    LocalFunctionImpPtrType lfptr_;

    // reference to local function 
    LocalFunctionImpType &lf_;

  public:
    //! Constructor initializing the underlying local function
    inline LocalFunctionWrapper ( const EntityType &entity,
                                  LocalFunctionStorageType &storage )
    : lfptr_( storage.getObject() ),
      lf_( *lfptr_ )
    {
      // init real local function with entity
      asImp().init( entity );
    }

    //! Constructor creating empty local function 
    inline explicit LocalFunctionWrapper ( LocalFunctionStorageType &storage ) 
    : lfptr_( storage.getObject() ),
      lf_( *lfptr_ )
    {
    }

    //! Constructor creating empty local function from given discrete
    //! function 
    template< class DiscreteFunctionType >
    inline explicit LocalFunctionWrapper ( DiscreteFunctionType &discreteFunction ) 
    : lfptr_( discreteFunction.localFunctionStorage().getObject() ),
      lf_( *lfptr_ )
    {
    }

    /** \brief copy constructor (does a shallow copy)
      
        \note this copy constructor must only called once, 
              i.e. on call of the method localFunction on the DiscreteFunction class,
              since this does only a shallow copy
     */
    inline LocalFunctionWrapper ( const ThisType &other )
    : lfptr_( other.lfptr_ ),
      lf_( *lfptr_ )
    {
      // make sure that we only got one reference of the lf_ object
      // otherwise it is not clear that changes will be correct 
      // if this assertion fails an unauthorized copy was done 
      // (i.e. copy-constructing a local function)
       
      assert( lfptr_.referenceCounter() == 1 ); 
    }

  private:
    // prohibit assignment 
    inline ThisType &operator= ( const ThisType & );

  protected:
    const LocalFunctionImpType &asImp () const
    { 
      return lf_;
    } 

    LocalFunctionImpType &asImp () 
    {
      return lf_;
    } 
  }; // end LocalFunctionWrapper  

  

  template< class LocalFunctionFactoryImp >
  class LocalFunctionStack
  : public Fem :: ObjectStack< LocalFunctionFactoryImp >
  {
  public:
    typedef LocalFunctionFactoryImp LocalFunctionFactoryType;

  private:
    typedef LocalFunctionStack< LocalFunctionFactoryType > ThisType;
    typedef Fem :: ObjectStack< LocalFunctionFactoryImp > BaseType;

  public:
    typedef LocalFunctionWrapper< ThisType > LocalFunctionType;

  public:
    inline explicit LocalFunctionStack ( const LocalFunctionFactoryType &factory )
    : BaseType( factory )
    {
    }

  private:
    inline LocalFunctionStack ( const ThisType & );

    inline ThisType &operator= ( const ThisType & );

  public:
    LocalFunctionType localFunction ()
    {
      return LocalFunctionType( *this );
    }
    
    template< class EntityType >
    const LocalFunctionType localFunction ( const EntityType &entity ) const
    {
      return LocalFunctionType( entity, *this );
    }

    template< class EntityType >
    LocalFunctionType localFunction ( const EntityType &entity )
    {
      return LocalFunctionType( entity, *this );
    }
  };
  
  } // end namespace Fem

} // end namespace Dune
#endif
