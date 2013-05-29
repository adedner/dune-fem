#ifndef DUNE_FEM_CODIMINDEXSET_HH
#define DUNE_FEM_CODIMINDEXSET_HH

#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>

#include <dune/fem/io/streams/xdrstreams.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/grid/utility/persistentcontainervector.hh>
#include <dune/grid/utility/persistentcontainerwrapper.hh>
#include <dune/grid/utility/persistentcontainermap.hh>

#ifdef ENABLE_ADAPTIVELEAFINDEXSET_FOR_YASPGRID
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

namespace Dune 
{

  // PersistentContainer for YaspGrid
  // --------------------------------

  template< int dim, class Data >
  class PersistentContainer< YaspGrid< dim >, Data >
  : public PersistentContainerVector< YaspGrid< dim >, 
                                      typename YaspGrid< dim >::LeafIndexSet,
                                      std::vector<Data> >
  {
    typedef YaspGrid< dim > Grid;
    typedef PersistentContainerVector< Grid, typename Grid::LeafIndexSet, std::vector<Data> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const Grid &grid, const int codim, const Data& value = Data() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.0, value )
    {}
  };

  // PersistentContainer for SGrid
  // -------------------------------

  template< int dim, int dimworld, class ctype, class Data >
  class PersistentContainer< SGrid< dim, dimworld, ctype >, Data >
  : public PersistentContainerVector< SGrid< dim, dimworld, ctype >, 
                                      typename SGrid< dim, dimworld, ctype >::LeafIndexSet,
                                      std::vector<Data> >
  {
    typedef SGrid< dim, dimworld, ctype > Grid ;
    typedef PersistentContainerVector< Grid, typename Grid::LeafIndexSet, std::vector<Data> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const Grid &grid, const int codim, const Data& value = Data() )
    : BaseType( grid, codim, grid.leafIndexSet(), 1.0, value )
    {}
  };

}
#endif


namespace Dune 
{

  namespace Fem
  {

    //***********************************************************************
    //
    //  Index Set for one codimension
    //  --CodimIndexSet 
    //
    //***********************************************************************
    template <class GridImp>  
    class CodimIndexSet
    {
    protected:  
      typedef GridImp GridType;
      typedef CodimIndexSet < GridType >  ThisType;

    private:
      enum INDEXSTATE { UNUSED = 0,  // unused indices
                        USED   = 1,  // used indices
                        NEW    = 2 };//  new indices 

      // reference to grid 
      const GridType& grid_;

    public:
      // type of exported index 
      typedef int IndexType ;

    protected:  
      // indices in this status have not been initialized 
      static IndexType invalidIndex() { return -1; }

      typedef std::pair< IndexType, INDEXSTATE > IndexPair;

      // array type for indices 
      typedef MutableArray< IndexType > IndexArrayType;

      class IndexPersistentContainer 
        : public PersistentContainer< GridImp, IndexPair >
      {
        typedef PersistentContainer< GridImp, IndexPair > BaseType ;

        template <class G, class T> 
        struct PublicPersistentContainerWrapper : public PersistentContainerWrapper< G, T >
        {
          using PersistentContainerWrapper< G, T > :: hostContainer_;
        };

        template <class G, class I, class V> 
        struct PublicPersistentContainerVector : public PersistentContainerVector< G, I, V >
        {
          using PersistentContainerVector< G, I, V > :: indexSet;
        };

        template <class G, class I, class M> 
        struct PublicPersistentContainerMap : public PersistentContainerMap< G, I, M >
        {
          using PersistentContainerMap< G, I, M > :: idSet;
        };

      public:
        using BaseType :: size ;
        using BaseType :: resize ;
        using BaseType :: codimension ;

        typedef typename BaseType :: Value Value;
        typedef typename BaseType :: Size  Size ;

        IndexPersistentContainer( const GridImp& grid, const int codim, const Value& value )
          : BaseType( grid, codim, value ) 
        {}

        void enlarge( const Value& value = Value() ) 
        {
          // call corrected implementation 
          enlargeImpl( *this, value ); 
        }

      protected:  
        template <class G, class T> 
        void enlargeImpl( PersistentContainerWrapper< G, T >& container, const Value& value ) 
        {
          enlargeImpl( ((PublicPersistentContainerWrapper< G, T > &) container).hostContainer_, value );
        }

        // enlarge implementation for persistent containers based on vectors 
        template < class G, class IndexSet, class Vector >
        void enlargeImpl( PersistentContainerVector< G, IndexSet, Vector >& container, const Value& value ) 
        {
          const Size indexSetSize = ((PublicPersistentContainerVector< G, IndexSet,
                Vector >& ) container).indexSet().size( codimension() );
          if( size() < indexSetSize ) 
            resize( value ); 
        }

        // enlarge implementation for persistent containers based on maps
        template < class G, class IdSet, class Map >
        void enlargeImpl( PersistentContainerMap< G, IdSet, Map >& container, const Value& value ) 
        {
          const Size idSetSize = ((PublicPersistentContainerMap< G, IdSet, Map >& ) container).idSet().size( codimension() );
          if( size() < idSetSize ) 
            resize( value ); 
        }
      };

      typedef IndexPersistentContainer IndexContainerType ;

      //typedef PersistentContainer< GridImp, IndexPair > IndexContainerType;

      // the mapping of the global to leaf index 
      IndexContainerType leafIndex_;

      // stack for holes 
      IndexArrayType holes_; 
     
      // Array that only remeber the occuring 
      // holes (for compress of data)
      IndexArrayType oldIdx_; 
      IndexArrayType newIdx_; 
     
      // next index to give away 
      IndexType nextFreeIndex_;

      // last size of set before compress (needed in parallel runs) 
      IndexType lastSize_;

      // codim for which index is provided 
      const int myCodim_; 

      // actual number of holes 
      IndexType numberHoles_;

    public:

      //! Constructor taking memory factor (default = 1.1)
      CodimIndexSet (const GridType& grid, 
                     const int codim, 
                     const double memoryFactor = 1.1) 
        : grid_( grid ) 
        , leafIndex_( grid, codim, IndexPair( invalidIndex(), UNUSED ) )
        , holes_(0)
        , oldIdx_(0)
        , newIdx_(0)
        , nextFreeIndex_ (0)
        , lastSize_ (0)
        , myCodim_( codim ) 
        , numberHoles_(0)
      {
        setMemoryFactor(memoryFactor);
      }

      //! set memory overestimation factor 
      void setMemoryFactor(const double memoryFactor)
      {
        holes_.setMemoryFactor(memoryFactor);
        oldIdx_.setMemoryFactor(memoryFactor);
        newIdx_.setMemoryFactor(memoryFactor);
      }

      //! reallocate the vectors
      void resize ()
      {
        // enlarge index container, do not shrink, because the old indices are still needed
        leafIndex_.enlarge( IndexPair( invalidIndex(), UNUSED ) );
      }

      //! prepare for setup (nothing to do here)
      void prepareCompress () {}

    public:  
      //! clear set 
      void clear() 
      {
        // set all values to invalidIndex  
        leafIndex_.fill( IndexPair( invalidIndex(), UNUSED ) );
        // reset next free index 
        nextFreeIndex_ = 0;
      }

      //! set all entries to unused 
      void resetUsed() 
      {
        typedef typename IndexContainerType::Iterator Iterator;
        const Iterator end = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != end; ++it )
          it->second = UNUSED;
      }

      bool consecutive ()
      {
        typedef typename IndexContainerType::Iterator Iterator;
        bool consecutive = true;
        const Iterator end = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != end; ++it )
          consecutive &= (it->first < nextFreeIndex_);
        return consecutive;
      }

      //! set all entries to unused 
      void checkConsecutive () { assert( consecutive() ); }

      //! clear holes, i.e. set number of holes to zero 
      void clearHoles() 
      {
        // set number of holes to zero 
        numberHoles_ = 0;
        // remember actual size 
        lastSize_ = nextFreeIndex_;
      }

      //! make to index numbers consecutive 
      //! return true, if at least one hole was closed 
      bool compress ()
      {
        const int sizeOfVecs = leafIndex_.size();
        holes_.resize( sizeOfVecs );

        // true if a least one dof must be copied 
        bool haveToCopy = false;
        
        // mark holes 
        int actHole = 0;
        int newActSize = 0;
        typedef typename IndexContainerType::Iterator Iterator;
        const Iterator end = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != end; ++it )
        {
          const IndexPair &leafIdx = *it;
          if( leafIdx.first != invalidIndex() )
          {
            // create vector with all holes 
            if( leafIdx.second == UNUSED )
            {
              holes_[actHole] = leafIdx.first;
              ++actHole;
            }

            // count the size of the leaf indices 
            ++newActSize;
          }
        }

        assert( newActSize >= actHole );
        // the new size is the actual size minus the holes 
        int actSize = newActSize - actHole;

        // resize hole storing vectors 
        oldIdx_.resize(actHole);
        newIdx_.resize(actHole);

        // only compress if number of holes > 0    
        if(actHole > 0)
        {
          // close holes 
          //
          // NOTE: here the holes closing should be done in 
          // the opposite way. future work. 
          int holes = 0; // number of real holes 
          //size_t i = 0;
          const Iterator end = leafIndex_.end();
          for( Iterator it = leafIndex_.begin(); it != end; ++it )
          {
            IndexPair &leafIdx = *it;
            // a index that is used but larger then actual size 
            // has to move to a hole 
            if( leafIdx.second == UNUSED ) 
            {
              // all unused indices are reset to invalidIndex  
              leafIdx.first = invalidIndex();
            }
            else 
            {
              // if used index lies behind size, then index has to move 
              // to one of the holes 
              if(leafIdx.first >= actSize)
              {
                // serach next hole that is smaler than actual size 
                actHole--;
                // if actHole < 0 then error, because we have index larger then
                // actual size 
                assert(actHole >= 0);
                while ( holes_[actHole] >= actSize )
                {
                  actHole--;
                  if(actHole < 0) break;
                }

                assert(actHole >= 0);

#if HAVE_MPI 
                // only for none-ghost elements hole storage is applied
                // this is because ghost indices might have in introduced 
                // after the resize was done. 
                if( leafIdx.second == USED ) 
#endif
                {
                  // remember old and new index 
                  oldIdx_[holes] = leafIdx.first; 
                  newIdx_[holes] = holes_[actHole];
                  ++holes;
                }
                
                leafIdx.first = holes_[actHole];

                // means that dof manager has to copy the mem
                leafIdx.second = NEW;
                haveToCopy = true;
              }
            }
          }

          // this call only sets the size of the vectors 
          oldIdx_.resize(holes);
          newIdx_.resize(holes);
        } // end if actHole > 0  
       
        // store number of actual holes 
        numberHoles_ = oldIdx_.size();

        // adjust size
        leafIndex_.resize( IndexPair( invalidIndex(), UNUSED ) );
        // shrinkToFit does not do anything 
        // leafIndex_.shrinkToFit();
        
        // the next index that can be given away is equal to size
        nextFreeIndex_ = actSize;

#ifndef NDEBUG
        checkConsecutive();
#endif

        return haveToCopy;
      }

      //! return how much extra memory is needed for restriction 
      IndexType additionalSizeEstimate () const { return nextFreeIndex_; }

      //! return size of grid entities per level and codim 
      IndexType size () const
      {
        return nextFreeIndex_;
      }
      
      //! return size of grid entities per level and codim 
      IndexType realSize () const
      {
        return leafIndex_.size();
      }

      //! return leaf index for given entity   
      //- --index 
      template <class EntityType>
      IndexType index ( const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        assert( checkValidIndex( leafIndex_[ entity ].first ) );
        return leafIndex_[ entity ].first;
      }
      
      //! return leaf index for given entity   
      template <class EntityType>
      IndexType subIndex ( const EntityType& entity,
                           const int subNumber ) const 
      {
        assert( 0 == EntityType :: codimension );
        assert( checkValidIndex( leafIndex_( entity, subNumber ).first ) );
        return leafIndex_( entity, subNumber ).first;
      }
      
      //! return state of index for given hierarchic number  
      template <class EntityType> 
      bool exists ( const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        return leafIndex_[ entity ].second != UNUSED;
      }
     
      template <class EntityType> 
      bool exists ( const EntityType& entity ,
                    const int subNumber ) const 
      {
        assert( 0 == EntityType :: codimension );
        return leafIndex_( entity, subNumber).second != UNUSED;
      }
     
      //! return number of holes 
      IndexType numberOfHoles () const
      {
        return numberHoles_;
      }

      //! return old index, for dof manager only 
      IndexType oldIndex (int elNum ) const
      {
        assert( numberHoles_ == IndexType(oldIdx_.size()) );
        return oldIdx_[elNum]; 
      }

      //! return new index, for dof manager only returns index 
      IndexType newIndex (int elNum) const
      {
        assert( numberHoles_ == IndexType(newIdx_.size()) );
        return newIdx_[elNum]; 
      }

      // insert element and create index for element number 
      template <class EntityType> 
      void insert (const EntityType& entity )
      {
        assert( myCodim_ == EntityType :: codimension );
        insertIdx( leafIndex_[ entity ] );
      }

      // insert element and create index for element number 
      template <class EntityType> 
      void insertSubEntity (const EntityType& entity, 
                            const int subNumber)  
      {
        assert( 0 == EntityType :: codimension );
        insertIdx( leafIndex_( entity, subNumber ) );
      }

      // insert element as ghost and create index for element number 
      template <class EntityType> 
      void insertGhost (const EntityType& entity )
      {
        assert( myCodim_ == EntityType :: codimension );
        // insert index 
        IndexPair &leafIdx = leafIndex_[ entity ];
        insertIdx( leafIdx );

        // if index is also larger than lastSize
        // mark as new to skip old-new index lists 
        if( leafIdx.first >= lastSize_ ) 
        {
          leafIdx.second = NEW;
        }
      }

      // insert element and create index for element number 
      template <class EntityType> 
      void markForRemoval( const EntityType& entity )
      {
        assert( myCodim_ == EntityType :: codimension );
        leafIndex_[ entity ].second = UNUSED;
      }

      // insert element as ghost and create index for element number 
      template <class EntityType> 
      bool validIndex (const EntityType& entity ) const
      {
        assert( myCodim_ == EntityType :: codimension );
        return leafIndex_[ entity ].first >= 0; 
      }

      void print( std::ostream& out ) const 
      {
        typedef typename IndexContainerType::ConstIterator Iterator;
        const Iterator end = leafIndex_.end();
        for( Iterator it = leafIndex_.begin(); it != end; ++it )
        {
          const IndexPair &leafIdx = *it;
          out << "idx: " << leafIdx.first << "  stat: " << leafIdx.second << std::endl;
        }
      }

    protected:
      // return true if the index idx is valid 
      bool checkValidIndex( const IndexType& idx ) const 
      {
        assert( idx != invalidIndex() );
        assert( idx  < size() );
        return (idx != invalidIndex() ) && ( idx < size() );
      }

      // insert element and create index for element number  
      void insertIdx ( IndexPair &leafIdx )
      {
        if( leafIdx.first == invalidIndex() )
          leafIdx.first = nextFreeIndex_++;
        leafIdx.second = USED;
      }

    public:  
      // write to stream 
      template <class StreamTraits> 
      bool write(OutStreamInterface< StreamTraits >& out) const
      {
        // store current index set size 
        out << nextFreeIndex_ ;
        
        // for consistency checking, write size as 64bit integer
        const uint64_t mysize = leafIndex_.size();
        out << mysize ;

        // backup indices 
        typedef typename IndexContainerType::ConstIterator ConstIterator;
        const ConstIterator end = leafIndex_.end();
        for( ConstIterator it = leafIndex_.begin(); it != end; ++it )
          out << (*it).first ;

        return true;
      }
      
      // read from stream 
      template <class StreamTraits> 
      bool read(InStreamInterface< StreamTraits >& in)
      {
        // read current index set size 
        in >> nextFreeIndex_ ;
        
        // for consistency checking 
        uint64_t storedSize = 0;
        in >> storedSize ;

        uint64_t leafsize = leafIndex_.size();
        // the stored size can be larger (visualization of parallel grids in serial)
        if( storedSize < leafsize ) 
        {
          DUNE_THROW(InvalidStateException,"CodimIndexSet: size consistency check failed during restore!"); 
        }

        // restore indices  
        typedef typename IndexContainerType::Iterator Iterator;
        const Iterator end = leafIndex_.end();
        uint64_t count = 0 ;
        for( Iterator it = leafIndex_.begin(); it != end; ++it, ++count )
          in >> (*it).first ;

        // also read indices that were stored but are not needed on read 
        if( count < storedSize )
        {
          IndexType value ;
          const uint64_t leftOver = storedSize - count ;
          for( uint64_t i = 0; i < leftOver; ++i ) 
            in >> value ;
        }

        return true;
      }
      
    }; // end of CodimIndexSet  

  } // namespace Fem

} // namespace Dune 

#endif // #ifndef DUNE_FEM_CODIMINDEXSET_HH
