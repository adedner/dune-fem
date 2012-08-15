#ifndef BASEFUNCTION_STORAGE_INTERFACE_HH
#define BASEFUNCTION_STORAGE_INTERFACE_HH

//- system includes 
#include <iostream>
#include <map>
#include <vector>
#include <list>

//- Dune includes 
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/fem/misc/threadmanager.hh>

namespace Dune 
{

  namespace Fem 
  {

    //! \brief Storage policy for base function sets.
    //! In a base function set, the base function values on quadrature points
    //! can either be cached or always recalculated. The storage policies
    //! CachingStorage and SimpleStorage do exactly that. The present class
    //! implements the common functionality and can be seen as a layer of
    //! abstraction in the access to basefunctions.
    template <int dimworld> 
    class StorageInterface
    {
        typedef StorageInterface<dimworld> ThisType;
        typedef std::list< ThisType* > StorageInterfaceListType;

        typedef StorageInterfaceListType* StorageInterfaceListPointer;
        typedef std::vector< size_t > QuadratureIdentifierType;
        typedef std::list< QuadratureIdentifierType > QuadratureListType;

        enum { ids = 0, codims = 1, sizes = 2, geoIndex = 3, sizeIndents = 4 };

        // return reference to list singleton pointer 
        static StorageInterfaceListPointer& storageListPtr () 
        {
          static StorageInterfaceListType* storageListObj = 0;
          return storageListObj;
        }

        // check if list is empty and if yes delete list 
        void checkAndDeleteStorageList() 
        {
          if( storageListPtr()->empty () )
          {
            delete storageListPtr();
            storageListPtr() = 0;
          }
        }

        // singelton implementation 
        static StorageInterfaceListType & storageList ()
        {
          // if list pointer is 0 then create new object 
          if( ! storageListPtr() )
          {
            // make sure we are in single thread mode 
            assert( ThreadManager :: singleThreadMode() );
            storageListPtr() = new StorageInterfaceListType();
          }

          assert( storageListPtr() );
          return *(storageListPtr());
        }

        // singelton implementation 
        static QuadratureListType & quadratureList ()
        {
          static QuadratureListType quadratureListObj;
          return quadratureListObj;
        }
        
      public:
        //! Constructor, add me to the list of storages 
        StorageInterface()
        {
          storageList().push_back(this);
        }

        //! Destructor, remove me from the list of storages 
        virtual ~StorageInterface() 
        {
          // make sure we are in single thread mode 
          assert( ThreadManager :: singleThreadMode() );
          typedef typename StorageInterfaceListType::iterator IteratorType;
          IteratorType endit = storageList().end();
          for(IteratorType it = storageList().begin(); it != endit; ++it)
          {
            if( (*it) == this )
            {
              storageList().erase(it);
              break;
            }
          }

          // if list is empty, then delete list 
          checkAndDeleteStorageList();
        }

        //! for a newly created storage cache all existing quadratures 
        template <class StorageImp>
        void cacheExistingQuadratures(StorageImp & storage)
        {
          // make sure we are in single thread mode 
          assert( ThreadManager :: singleThreadMode() );

          // get geometry index of storage  
          const unsigned int storageIndex = GlobalGeometryTypeIndex :: index( storage.geometryType() );

          typedef typename QuadratureListType::iterator IteratorType;
          IteratorType endit = quadratureList().end();
          for(IteratorType it = quadratureList().begin(); it != endit; ++it)
          {
            // only cache base functions for quadratures with same geometry type 
            if ( (*it)[geoIndex] == storageIndex )
            {
              // get if and codim of quad 
              const size_t id = (*it)[ ids ];
              const size_t codim = (*it)[ codims ];
              const size_t quadSize = (*it)[ sizes ];
              storage.cacheQuadrature(id, codim, quadSize);
            }
          }
        }

        //! cache quadrature for given id and codim 
        virtual void cacheQuadrature(const size_t id, 
                                     const size_t codim, 
                                     const size_t quadSize ) const = 0;

        //! return geometry type of base function set 
        virtual GeometryType geometryType() const = 0;

        template <class QuadratureType>
        static void registerQuadratureToStorages(const QuadratureType & quad)
        {
          const size_t codim = QuadratureType :: codimension;
          registerQuadratureToStorages(quad, quad.geometryType(), codim);
        }

        //! register quadrature for all existing storages 
        template <class QuadratureType>
        static void registerQuadratureToStorages(const QuadratureType & quad, 
                                                 const GeometryType& geoType,   
                                                 const size_t codim)
        {
          // make sure we are in single thread mode 
          assert( ThreadManager :: singleThreadMode() );

          const size_t id = quad.id();
          const size_t quadSize = quad.nop();
          // store quadrature 
          QuadratureIdentifierType ident( sizeIndents );

          ident[ ids ]    = id ;
          ident[ codims ] = codim; 
          ident[ sizes ]  = quadSize;
          ident[ geoIndex ] = GlobalGeometryTypeIndex :: index( geoType ) ;
          quadratureList().push_back(ident);

          typedef typename StorageInterfaceListType::iterator IteratorType;
          IteratorType endit = storageList().end();
          for(IteratorType it = storageList().begin(); it != endit; ++it)
          {
            // make sure that only base functions and quadratures 
            // with same geometry type are cached 
            if ( geoType == (*it)->geometryType() )
              (*it)->cacheQuadrature(id, codim, quadSize);
          }
        }
    };

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef BASEFUNCTION_STORAGE_INTERFACE_HH
