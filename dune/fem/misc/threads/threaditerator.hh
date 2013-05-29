#ifndef DUNE_FEM_THREADITERATOR_HH
#define DUNE_FEM_THREADITERATOR_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/misc/threadmanager.hh>

namespace Dune 
{

  namespace Fem 
  {

    /** \brief Thread iterator */
    template <class DiscreteFunctionSpace>  
    class ThreadIterator
    {
      typedef DiscreteFunctionSpace SpaceType;
    public:  
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename SpaceType :: IteratorType IteratorType ;
      typedef typename IteratorType :: Entity EntityType ;
      typedef typename SpaceType :: GridType :: template Codim<0> ::
        EntityPointer EntityPointer ;
      typedef typename SpaceType :: IndexSetType IndexSetType ;
    protected:  
      const SpaceType& space_ ;
      const IndexSetType& indexSet_;

#ifdef USE_SMP_PARALLEL
      int sequence_;
      std::vector< IteratorType > iterators_;
      MutableArray< int > threadNum_;
      std::vector< std::vector< int > > threadId_; 
#endif
    public:  
      //! contructor creating thread iterators 
      explicit ThreadIterator( const SpaceType& spc )
        : space_( spc )
        , indexSet_( space_.indexSet() )
#ifdef USE_SMP_PARALLEL
        , sequence_( -1 )  
        , iterators_( ThreadManager::maxThreads() + 1 , space_.end() )
        , threadId_( ThreadManager::maxThreads() )
#endif
      {
#ifdef USE_SMP_PARALLEL
        threadNum_.setMemoryFactor( 1.1 ); 
#endif
        update();
      }

      //! return reference to space 
      const SpaceType& space() const { return space_; }

      //! update internal list of iterators 
      void update() 
      {
#ifdef USE_SMP_PARALLEL
        // if grid got updated also update iterators 
        if( sequence_ != space_.sequence() )
        {
          if( ! ThreadManager :: singleThreadMode() ) 
          {
            std::cerr << "Don't call ThreadIterator::update in a parallel environment!" << std::endl;
            assert( false );
            abort();
          }

          const size_t maxThreads = ThreadManager :: maxThreads() ;

          // get end iterator
          const IteratorType endit = space_.end();

          IteratorType it = space_.begin(); 
          if( it == endit ) 
          {
            // set all iterators to end iterators 
            for( size_t thread = 0; thread <= maxThreads; ++thread ) 
              iterators_[ thread ] = endit ;

            // free memory here 
            threadNum_.resize( 0 );

            // update sequence number 
            sequence_ = space_.sequence();
            return ;
          }

          // thread 0 starts at begin 
          iterators_[ 0 ] = it ;

          // get size for index set 
          const size_t size = indexSet_.size( 0 );

          // resize threads storage 
          threadNum_.resize( size );
          // set all values to default value 
          for(size_t i = 0; i<size; ++i) threadNum_[ i ] = -1;

          // here use iterator to count 
          size_t checkSize = 0;
          const size_t roundOff = (size % maxThreads);
          const size_t counterBase = ((size_t) size / maxThreads );
          for( size_t thread = 1; thread <= maxThreads; ++thread ) 
          {
            size_t i = 0; 
            const size_t counter = counterBase + (( (thread-1) < roundOff ) ? 1 : 0);
            checkSize += counter ;
            //std::cout << counter << " for thread " << thread-1 << std::endl;
            while( (i < counter) && (it != endit) )
            {
              assert( indexSet_.index( *it ) < (size_t) threadNum_.size() );
              threadNum_[ indexSet_.index( *it ) ] = thread - 1;
              ++i;
              ++it;
            }
            iterators_[ thread ] = it ;
          }
          iterators_[ maxThreads ] = endit ;

          if( checkSize != size ) 
          {
            assert( checkSize == size );
            DUNE_THROW(InvalidStateException,"Partitioning inconsistent!"); 
          }

          // update sequence number 
          sequence_ = space_.sequence();

          //for(size_t i = 0; i<size; ++i ) 
          //  std::cout << threadNum_[ i ] << std::endl;
        }
#endif
      }

      //! return begin iterator for current thread 
      IteratorType begin() const 
      {
#ifdef USE_SMP_PARALLEL
        return iterators_[ ThreadManager :: thread() ];
#else 
        return space_.begin();
#endif
      }

      //! return end iterator for current thread 
      IteratorType end() const 
      {
#ifdef USE_SMP_PARALLEL
        return iterators_[ ThreadManager :: thread() + 1 ];
#else 
        return space_.end();
#endif
      }

      //! return thread number this entity belongs to 
      int index(const EntityType& entity ) const 
      {
        return indexSet_.index( entity );
      }

      //! return thread number this entity belongs to 
      int thread(const EntityType& entity ) const 
      {
#ifdef USE_SMP_PARALLEL
        assert( (size_t) threadNum_.size() > indexSet_.index( entity ) );
        // NOTE: this number can also be negative for ghost elements or elements
        // that do not belong to the set covered by the space iterators 
        return threadNum_[ indexSet_.index( entity ) ];
#else 
        return 0;
#endif
      }
    };


    /** \brief Thread iterator */
    template <class DiscreteFunctionSpace>  
    class ThreadIteratorPointer
    {
      typedef DiscreteFunctionSpace SpaceType;
    public:  
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;
      typedef typename SpaceType :: GridType :: template Codim<0> :: Entity EntityType ;
      typedef typename SpaceType :: GridType :: template Codim<0> ::
        EntityPointer EntityPointer ;
      typedef typename SpaceType :: IndexSetType IndexSetType ;
    protected:  
      const SpaceType& space_ ;
      const IndexSetType& indexSet_;

      class Iterator
      {
        const std::vector< EntityPointer >& vec_;
        size_t pos_ ;
      public:
        typedef typename SpaceType :: GridType :: template Codim<0> :: Entity EntityType ;
        typedef EntityType Entity ;
        Iterator( const std::vector< EntityPointer >& vec, const size_t pos)
          : vec_( vec ), pos_( pos )
        {}

        const EntityType& operator* () const 
        {
          assert( pos_ < vec_.size() );
          return *( vec_[ pos_ ] );
        }

        void operator ++ () { ++pos_; }
        bool operator != (const Iterator & other ) const 
        {
          return pos_ != other.pos_ ;
        }
      };

      int sequence_;
      std::vector< std::vector< EntityPointer > > pointers_;
      std::vector< int > threadNum_;
    public:  
      typedef Iterator IteratorType ;
      explicit ThreadIteratorPointer( const SpaceType& spc )
        : space_( spc )
        , indexSet_( space_.indexSet() )
#ifdef USE_SMP_PARALLEL
        , sequence_( -1 )  
        , pointers_( ThreadManager :: maxThreads() + 1 )
#endif
      {
        update();
      }

      void update() 
      {
#ifdef USE_SMP_PARALLEL
        if( sequence_ != space_.sequence() )
        {
          if( ! ThreadManager :: singleThreadMode() ) 
          {
            std::cerr << "Don't call ThreadIterator::update in a parallel environment!" << std::endl;
            assert( false );
            abort();
          }

          const int maxThreads = ThreadManager :: maxThreads() ;
          typedef typename SpaceType :: IteratorType SpcIteratorType ;
          SpcIteratorType it = space_.begin(); 
          const SpcIteratorType endit = space_.end();
          if( it == endit ) 
          {
            // update sequence number 
            sequence_ = space_.sequence();
            return ;
          }

          int size = 0;
          for( SpcIteratorType countit = space_.begin(); countit != endit; ++countit, ++size ) {}
          // resize threads 
          threadNum_.resize( size );

          // here ruse iterator to count 
          const int counter = ((int) size / maxThreads );
          for( int thread = 0; thread <= maxThreads; ++thread ) 
          {
            int i = 0; 
            pointers_[ thread ].resize( counter, it );
            while( (i < counter) && (it != endit) )
            {
              pointers_[ thread ][ i ] = it ;
              threadNum_[ indexSet_.index( *it ) ] = thread ;
              ++i;
              ++it;
            }
          }
          // update sequence number 
          sequence_ = space_.sequence();
        }
#endif
      }

      Iterator begin() const 
      {
#ifdef USE_SMP_PARALLEL
        return Iterator( pointers_[ ThreadManager::thread() ], 0 );
#else 
        return Iterator( pointers_[ 0 ], 0 );
#endif
      }

      Iterator end() const 
      {
#ifdef USE_SMP_PARALLEL
        return Iterator( pointers_[ ThreadManager::thread() ], 
                         pointers_[ ThreadManager::thread() ].size() );
#else 
        return Iterator( pointers_[ 0 ], pointers_[ 0 ].size() );
#endif
      }

      int thread(const EntityType& entity ) const 
      {
        assert( ((int) threadNum_.size()) > indexSet_.index( entity ) );
        return threadNum_[ indexSet_.index( entity ) ];
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_THREADITERATOR_HH
