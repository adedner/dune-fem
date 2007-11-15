#ifndef DUNE_DISCRETEFUNCTION_CC
#define DUNE_DISCRETEFUNCTION_CC

#include <fstream>

#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/streams/streams.hh>

namespace Dune 
{
 
  // Default Implementations 
  // -----------------------

  template< class DiscreteFunctionTraits >
  inline void DiscreteFunctionDefault< DiscreteFunctionTraits > :: clear ()
  {
    const DofIteratorType end = dend();
    for( DofIteratorType it = dbegin(); it != end; ++it )
      *it = 0;
  }

  template< class Traits >
  inline void DiscreteFunctionDefault< Traits >
    :: addScaled ( const DiscreteFunctionType &g,
                   const RangeFieldType &s )
  {
    assert( this->size() == g.size() );
    const DofIteratorType end = this->dend();
    ConstDofIteratorType git = g.dbegin();
    for( DofIteratorType it = this->dbegin(); it != end; ++it, ++git )
      (*it) += s * (*git);
  }

  template< class Traits >
  inline typename DiscreteFunctionDefault< Traits > :: RangeFieldType *
  DiscreteFunctionDefault< Traits > :: allocDofPointer()
  {
    dofPointerLock_.lock();

    const unsigned int size = this->size();
    RangeFieldType *dofPointer = new RangeFieldType[ size ];
    
    unsigned int i = 0;
    const DofIteratorType end = dend();
    for( DofIteratorType it = dbegin(); it != end; ++it )
      dofPointer[ i++ ] = *it;
    assert( i == size );

    return dofPointer;
  }

  template< class Traits >
  inline void DiscreteFunctionDefault< Traits >
    :: freeDofPointer( RangeFieldType *dofPointer )
  {
    unsigned int i = 0;
    const DofIteratorType end = dend();
    for( DofIteratorType it = dbegin(); it != end; ++it )
      *it = dofPointer[ i++ ];
    assert( i == size() );

    delete[] dofPointer;
    dofPointerLock_.unlock();
  }

// scalarProductDofs
template <class DiscreteFunctionTraits>
inline typename DiscreteFunctionTraits::DiscreteFunctionSpaceType::RangeFieldType 
DiscreteFunctionDefault<DiscreteFunctionTraits>::
scalarProductDofs(const DiscreteFunctionType& g) const
{
  typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType; 
  assert(this->size() == g.size());

  RangeFieldType skp = 0.;

  ConstDofIteratorType endit = this->dend ();
  ConstDofIteratorType git =  g.dbegin ();

  // multiply
  for(ConstDofIteratorType it = this->dbegin(); it != endit; ++it,++git)
  {   
    skp += (*it) * (*git);
  }
  
  return skp;
}



  template< class DiscreteFunctionTraits >
  inline void
  DiscreteFunctionDefault<DiscreteFunctionTraits >
    :: assign( const DiscreteFunctionType &g )
  {
    assert( size() == g.size() );

    const DofIteratorType end = dend();
    ConstDofIteratorType git = g.dbegin();
    for( DofIteratorType it = dbegin(); it != end; ++it, ++git )
      *it = *git;
  }



// operator +=
/** \todo This operator can add a discretefunction defined on all levels to another
 * one defined only on one level.  We should somehow issue a warning in this case.
 */
template<class DiscreteFunctionTraits>
inline typename DiscreteFunctionDefault<DiscreteFunctionTraits> :: DiscreteFunctionType&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator += ( const DiscreteFunctionType& g ) 
{
  assert(this->size() == g.size());

  DofIteratorType endit = this->dend ();
  ConstDofIteratorType git = g.dbegin ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) 
  {
    *it += *git;
  }
  return asImp();
}

// operator -=
template<class DiscreteFunctionTraits>
inline typename DiscreteFunctionDefault<DiscreteFunctionTraits> :: DiscreteFunctionType&
DiscreteFunctionDefault<DiscreteFunctionTraits >::
operator -= ( const DiscreteFunctionType& g ) 
{
  assert(this->size() == g.size());

  DofIteratorType endit = this->dend ();
  ConstDofIteratorType git = g.dbegin ();
  for(DofIteratorType it = this->dbegin(); it != endit; ++it, ++git) 
  {
    *it -= *git;
  }
  return asImp();
}



  // operator *=
  template< class DiscreteFunctionTraits >
  inline
  typename DiscreteFunctionDefault< DiscreteFunctionTraits >
    :: DiscreteFunctionType &
  DiscreteFunctionDefault< DiscreteFunctionTraits >
    :: operator*= ( const RangeFieldType &scalar )
  {
    const DofIteratorType end = dend();
    for( DofIteratorType it = dbegin(); it != end; ++it )
      *it *= scalar;
    return asImp();
  }



  template< class DiscreteFunctionTraits >
  inline bool DiscreteFunctionDefault< DiscreteFunctionTraits >
    :: operator== ( const DiscreteFunctionType &g ) const
  {
    if( size() != g.size() )
      return false;
    
    const ConstDofIteratorType end = dend();

    ConstDofIteratorType fit = dbegin();
    ConstDofIteratorType git = g.dbegin();
    for( ; fit != end; ++fit, ++git )
      if( *fit != *git )
        return false;
    
    return true;
  }


  // print 
  template< class DiscreteFunctionTraits >
  inline void DiscreteFunctionDefault<DiscreteFunctionTraits >
    :: print ( std::ostream &out ) const
  {
    out << name() << std::endl;
    
    const ConstDofIteratorType end = dend();
    for( ConstDofIteratorType dit = dbegin(); dit != end; ++dit )
      out << (*dit) << std::endl;
  }



  template< class DiscreteFunctionTraits >
  template< class StreamTraits >
  inline void DiscreteFunctionDefault< DiscreteFunctionTraits >
    :: read ( InStreamInterface< StreamTraits > &in )
  {
    int sz;
    in >> sz;
    if( sz != size() )
      DUNE_THROW( IOError, "Trying to read discrete function of different size." );

    const DofIteratorType end = dend();
    for( DofIteratorType it = dbegin(); it != end; ++it )
      in >> *it;
  }



  template< class DiscreteFunctionTraits >
  template< class StreamTraits >
  inline void DiscreteFunctionDefault< DiscreteFunctionTraits >
    :: write ( OutStreamInterface< StreamTraits > &out ) const
  {
    out << size();

    const ConstDofIteratorType end = dend();
    for( ConstDofIteratorType it = dbegin(); it != end; ++it )
      out << *it;
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
  template< class DiscreteFunctionTraits >
  inline std :: ostream &
    operator<< ( std :: ostream &out,
                 const DiscreteFunctionInterface< DiscreteFunctionTraits > &df )
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
  template< class StreamTraits, class DiscreteFunctionTraits >
  inline OutStreamInterface< StreamTraits > &
    operator<< ( OutStreamInterface< StreamTraits > &out,
                 const DiscreteFunctionInterface< DiscreteFunctionTraits > &df )
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
   *  \returns the output stream (for concatenation)
   */
  template< class StreamTraits, class DiscreteFunctionTraits >
  inline InStreamInterface< StreamTraits > &
    operator>> ( InStreamInterface< StreamTraits > &in,
                 DiscreteFunctionInterface< DiscreteFunctionTraits > &df )
  {
    df.read( in );
    return in;
  }

} // end namespace Dune
#endif
