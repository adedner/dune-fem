#ifndef DUNE_FEM_DISCRETEFUNCTION_INLINE_HH
#define DUNE_FEM_DISCRETEFUNCTION_INLINE_HH

#include <fstream>

#include <dune/geometry/referenceelements.hh>

#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/misc/threadmanager.hh>
#include <dune/fem/gridpart/dunefemindexsets.hh>

#include "discretefunction.hh"

namespace Dune 
{
 
  namespace Fem 
  {

    // DiscreteFunctionDefault 
    // -----------------------

    template< class Traits >
    inline DiscreteFunctionDefault< Traits >
      :: DiscreteFunctionDefault ( const std::string &name,
                                   const DiscreteFunctionSpaceType &dfSpace,
                                   const LocalFunctionFactoryType &lfFactory )
    : dfSpace_( dfSpace ),
#ifdef USE_SMP_PARALLEL
      lfStorageVec_ ( Fem :: ThreadManager :: maxThreads() ),
#else 
      lfStorage_( lfFactory ),
#endif
      name_( name ),
      scalarProduct_( dfSpace )
    {
#ifdef USE_SMP_PARALLEL
      for( size_t i=0 ; i<lfStorageVec_.size(); ++i ) 
      {
        lfStorageVec_[ i ] = new LocalFunctionStorageType( lfFactory );
      }
#endif
    }

    template< class Traits >
    inline DiscreteFunctionDefault< Traits >::~DiscreteFunctionDefault ()
    {
      assert( !dofPointerLock_ );
#ifdef USE_SMP_PARALLEL
      for( size_t i=0 ; i<lfStorageVec_.size(); ++i )
      {
        delete lfStorageVec_[ i ]; 
      }
#endif
    }


    template< class Traits >
    inline const typename DiscreteFunctionDefault< Traits > ::  LocalFunctionType
    DiscreteFunctionDefault< Traits >
      :: localFunction ( const EntityType &entity ) const
    {
      return localFunctionStorage().localFunction( entity );
    }


    template< class Traits >
    inline typename DiscreteFunctionDefault< Traits > ::  LocalFunctionType
    DiscreteFunctionDefault< Traits >
      :: localFunction ( const EntityType &entity )
    {
      return localFunctionStorage().localFunction( entity );
    }


    template< class Traits >
    inline void DiscreteFunctionDefault< Traits > :: clear ()
    {
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        *it = 0;
    }

    
    template< class Traits >
    inline typename DiscreteFunctionDefault< Traits > :: RangeFieldType *
    DiscreteFunctionDefault< Traits > :: allocDofPointer()
    {
      dofPointerLock_.lock();

      const unsigned int size = BaseType :: size();
      RangeFieldType *dofPointer = new RangeFieldType[ size ];
      
      unsigned int i = 0;
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        dofPointer[ i++ ] = *it;
      assert( i == size );

      return dofPointer;
    }

    
    template< class Traits >
    inline void DiscreteFunctionDefault< Traits >
      :: freeDofPointer( RangeFieldType *dofPointer )
    {
      unsigned int i = 0;
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        *it = dofPointer[ i++ ];
      assert( i == BaseType :: size() );

      delete[] dofPointer;
      dofPointerLock_.unlock();
    }


    template< class Traits >
    inline void DiscreteFunctionDefault< Traits >
      ::axpy ( const RangeFieldType &s, const DiscreteFunctionInterfaceType &g )
    {
      assert( BaseType::size() == g.size() );
      const DofIteratorType end = BaseType::dend();
      ConstDofIteratorType git = g.dbegin();
      for( DofIteratorType it = BaseType::dbegin(); it != end; ++it, ++git )
        (*it) += s * (*git);
    }


    template< class Traits >
    inline typename DiscreteFunctionDefault< Traits > :: RangeFieldType
    DiscreteFunctionDefault< Traits >
      ::scalarProductDofs ( const DiscreteFunctionInterfaceType &other ) const
    {
      return scalarProduct_.scalarProductDofs( *this, other );
    }

    
    template< class Traits >
    inline void DiscreteFunctionDefault<Traits >
      :: print ( std::ostream &out ) const
    {
      out << BaseType :: name() << std::endl;
      
      const ConstDofIteratorType end = BaseType :: dend();
      for( ConstDofIteratorType dit = BaseType :: dbegin(); dit != end; ++dit )
        out << (*dit) << std::endl;
    }


    template< class Traits >
    inline bool DiscreteFunctionDefault< Traits >
      :: dofsValid () const
    {
      const ConstDofIteratorType end = BaseType :: dend();
      for( ConstDofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        if( *it != *it )
          return false;

      return true;
    }

    
    template< class Traits >
    inline void DiscreteFunctionDefault< Traits >
      ::assign ( const DiscreteFunctionInterfaceType &g )
    {
      assert( BaseType::size() == g.size() );

      const DofIteratorType end = BaseType::dend();
      ConstDofIteratorType git = g.dbegin();
      for( DofIteratorType it = BaseType::dbegin(); it != end; ++it, ++git )
        *it = *git;
    }


    template< class Traits >
    template< class Operation >
    inline typename DiscreteFunctionDefault< Traits >
      :: template CommDataHandle< Operation > :: Type
    DiscreteFunctionDefault< Traits > :: dataHandle ( const Operation *operation )
    {
      return BaseType :: space().createDataHandle( asImp(), operation );
    }


    template< class Traits >
    inline void DiscreteFunctionDefault< Traits >
      :: evaluate ( const DomainType &x,
                    RangeType &ret ) const
    {
      FieldVector< int, 0 > diffVariable;
      BaseType :: evaluate( diffVariable, x, ret );
    }


    template< class Traits >
    template< int diffOrder >
    inline void DiscreteFunctionDefault< Traits >
      ::evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                   const DomainType &x,
                   RangeType &ret ) const
    {
      typedef typename DiscreteFunctionSpaceType::IteratorType Iterator;
      typedef typename DiscreteFunctionSpaceType::EntityType   Entity;
      typedef typename Entity::Geometry Geometry;
      typedef typename Geometry :: LocalCoordinate LocalCoordinateType;

      const int dimLocal = LocalCoordinateType :: dimension;
      
      const DiscreteFunctionSpaceType &space = BaseType::space();
      const Iterator end = space.end();
      for( Iterator it = space.begin(); it != end; ++it )
      {
        const Entity &entity = *it;
        const Geometry &geometry = entity.geometry();

        const GenericReferenceElement< DomainFieldType, dimLocal > &refElement
          = GenericReferenceElements< DomainFieldType, dimLocal >::general( geometry.type() );

        const LocalCoordinateType xlocal = geometry.local( x );
        if( refElement.checkInside( xlocal ) )
        {
          const LocalFunctionType localFunction = BaseType::localFunction( entity );
          localFunction.evaluate( diffVariable, xlocal, ret );
          return;
        }
      }
      DUNE_THROW( RangeError, "DiscreteFunctionDefault::evaluate: x is not within domain." );
    }


    template< class Traits >
    inline typename DiscreteFunctionDefault< Traits > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Traits >
      ::operator+= ( const DiscreteFunctionInterfaceType &g )
    {
      assert( BaseType::size() == g.size() );

      const DofIteratorType end = BaseType::dend();
      ConstDofIteratorType git = g.dbegin();
      for( DofIteratorType it = BaseType::dbegin(); it != end; ++it, ++git )
        *it += *git;
      return asImp();
    }


    template< class Traits >
    template< class DFType >
    inline typename DiscreteFunctionDefault< Traits > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Traits >
      :: operator-= ( const DFType &g )
    {
      assert( BaseType :: size() == g.size() );

      const DofIteratorType end = BaseType :: dend();
      typename DFType :: ConstDofIteratorType git = g.dbegin();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it, ++git )
        *it -= *git;
      return asImp();
    }


    template< class Traits >
    inline typename DiscreteFunctionDefault< Traits > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Traits >
      :: operator*= ( const RangeFieldType &scalar )
    {
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        *it *= scalar;
      return asImp();
    }


    template< class Traits >
    inline typename DiscreteFunctionDefault< Traits > :: DiscreteFunctionType &
    DiscreteFunctionDefault< Traits >
      :: operator/= ( const RangeFieldType &scalar )
    {
      return BaseType :: operator*=( RangeFieldType( 1 ) / scalar );
    }


    template< class Traits >
    template< class StreamTraits >
    inline void DiscreteFunctionDefault< Traits >
      :: read ( InStreamInterface< StreamTraits > &in )
    {
      unsigned int versionId = in.readUnsignedInt();
      if( versionId < DUNE_VERSION_ID(0,9,1) )
        DUNE_THROW( IOError, "Trying to read outdated file." );
      else if( versionId > DUNE_MODULE_VERSION_ID(DUNE_FEM) )
        std :: cerr << "Warning: Reading discrete function from newer version: "
                    << versionId << std :: endl;

      // read name 
      in >> name_;

      // read size as integer 
      int mysize;
      in >> mysize ;
      
      // check size 
      if( mysize != BaseType :: size() && 
          BaseType :: size() != this->space().size() ) // only read compressed vectors 
      {
        DUNE_THROW( IOError, "Trying to read discrete function of different size." );
      }

      // read all dofs 
      const DofIteratorType end = BaseType :: dend();
      for( DofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        in >> *it;
    }


    template< class Traits >
    template< class StreamTraits >
    inline void DiscreteFunctionDefault< Traits >
      :: write ( OutStreamInterface< StreamTraits > &out ) const
    {
      unsigned int versionId = DUNE_MODULE_VERSION_ID(DUNE_FEM);
      out << versionId ;
      out << name_;
    
      // only allow write when vector is compressed 
      if( BaseType :: size() != this->space().size() )
        DUNE_THROW(InvalidStateException,"Writing DiscreteFunction in uncompressed state!");
      
      // write size as integer 
      const int mysize = BaseType :: size();
      out << mysize;

      // write all dofs 
      const ConstDofIteratorType end = BaseType :: dend();
      for( ConstDofIteratorType it = BaseType :: dbegin(); it != end; ++it )
        out << *it;
    }

    template< class Traits >
    bool DiscreteFunctionDefault< Traits >
      :: read_xdr ( const std :: string filename )
    {
      try
      {
        XDRFileInStream in( filename );
        BaseType :: read( in );
        return true;
      }
      catch( Exception e )
      {
        return false;
      }
    }
   
    
    template< class Traits >
    bool DiscreteFunctionDefault< Traits >
      :: write_xdr ( const std :: string filename ) const
    {
      try
      {
        XDRFileOutStream out( filename );
        BaseType :: write( out );
        return true;
      }
      catch( Exception e )
      {
        return false;
      }
    }
    
    
    template< class Traits >
    bool DiscreteFunctionDefault< Traits >
      :: read_ascii ( const std :: string filename )
    {
      try
      {
        ASCIIInStream in( filename );
        BaseType :: read( in );
        return true;
      }
      catch( Exception e )
      {
        return false;
      }
    }
   
    
    template< class Traits >
    bool DiscreteFunctionDefault< Traits >
      :: write_ascii ( const std :: string filename ) const
    {
      try
      {
        ASCIIOutStream out( filename );
        BaseType :: write( out );
        return true;
      }
      catch( Exception e )
      {
        return false;
      }
    }


    template< class Traits >
    void DiscreteFunctionDefault< Traits >
      :: backup() const
    {
      // get backup stream from persistence manager and write to it 
      write( PersistenceManager :: backupStream() );
    }


    template< class Traits >
    void DiscreteFunctionDefault< Traits >
      :: restore()
    {
      // get restore stream from persistence manager and read from it 
      read( PersistenceManager :: restoreStream() );
    }

    template< class Traits >
    void DiscreteFunctionDefault< Traits >
      :: insertSubData()
    {
      // if indexset is persistent it must be  
      // derived from PersistentIndexSetInterface 
      if( space().indexSet().persistent() )
      {
        // this marks the index set in the DofManager's list of index set as persistent
        PersistentIndexSetInterface& indexSet = (PersistentIndexSetInterface &) space().indexSet();
        indexSet.addBackupRestore();
      }
    }

    template< class Traits >
    void DiscreteFunctionDefault< Traits >
      :: removeSubData()
    {
      // if indexset is persistent it must be 
      // derived from PersistentIndexSetInterface 
      if( space().indexSet().persistent() )
      {
        // this unmarks the index set in the DofManager's list of index set as persistent
        PersistentIndexSetInterface& indexSet = (PersistentIndexSetInterface &) space().indexSet();
        indexSet.removeBackupRestore();
      }
    }

    template< class Traits >
    inline void DiscreteFunctionDefault< Traits >
      :: enableDofCompression ()
    {}


    template< class Traits >
    inline bool DiscreteFunctionDefault< Traits >
      :: operator== ( const DiscreteFunctionType &g ) const
    {
      if( BaseType :: size() != g.size() )
        return false;
      
      const ConstDofIteratorType end = BaseType :: dend();

      ConstDofIteratorType fit = BaseType :: dbegin();
      ConstDofIteratorType git = g.dbegin();
      for( ; fit != end; ++fit, ++git )
        if( *fit != *git )
          return false;
      
      return true;
    }
    
   
    template< class Traits >
    inline bool DiscreteFunctionDefault< Traits >
      :: operator!= ( const DiscreteFunctionType &g ) const
    {
      return !(operator==( g ));
    }

    
    template< class Traits >
    inline typename DiscreteFunctionDefault< Traits > :: LocalFunctionStorageType &
    DiscreteFunctionDefault< Traits > :: localFunctionStorage () const
    {
#ifdef USE_SMP_PARALLEL
      return *(lfStorageVec_[ Fem :: ThreadManager :: thread() ]);
#else 
      return lfStorage_;
#endif
    }



    // Stream Operators
    // ----------------

    /** \brief write a discrete function into an STL stream
     *  \relates DiscreteFunctionInterface
     *
     *  \param[in]  out  STL stream to write to
     *  \param[in]  df   discrete function to write
     *
     *  \returns the STL stream (for concatenation)
     */
    template< class Traits >
    inline std :: ostream &
      operator<< ( std :: ostream &out,
                   const DiscreteFunctionInterface< Traits > &df )
    {
      df.print( out );
      return out;
    }



    /** \brief write a discrete function into an output stream
     *  \relates DiscreteFunctionInterface
     *  \relatesalso OutStreamInterface
     *
     *  \param[in]  out  stream to write to
     *  \param[in]  df   discrete function to write
     *
     *  \returns the output stream (for concatenation)
     */
    template< class StreamTraits, class Traits >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out,
                   const DiscreteFunctionInterface< Traits > &df )
    {
      df.write( out );
      return out;
    }



    /** \brief read a discrete function from an input stream
     *  \relates DiscreteFunctionInterface
     *  \relatesalso InStreamInterface
     *
     *  \param[in]   in  stream to read from
     *  \param[out]  df  discrete function to read
     *
     *  \returns the input stream (for concatenation)
     */
    template< class StreamTraits, class Traits >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in,
                   DiscreteFunctionInterface< Traits > &df )
    {
      df.read( in );
      return in;
    }

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_DISCRETEFUNCTION_INLINE_HH
