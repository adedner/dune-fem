#ifndef DUNE_FEM_STREAMS_INLINE_HH
#define DUNE_FEM_STREAMS_INLINE_HH

#include <vector>

#include "streams.hh"

namespace Dune
{

  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const double value )
  {
    out.writeDouble( value );
    return out;
  }
  
  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const float value )
  {
    out.writeFloat( value );
    return out;
  }

  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const int value )
  {
    out.writeInt( value );
    return out;
  }

  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const char value )
  {
    out.writeChar( value );
    return out;
  }

  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const bool value )
  {
    out.writeBool( value );
    return out;
  }

  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const std :: string &s )
  {
    out.writeString( s );
    return out;
  }
  
  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const unsigned int value )
  {
    out.writeUnsignedInt( value );
    return out;
  }

  template< class Traits >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const unsigned long value )
  {
    out.writeUnsignedLong( value );
    return out;
  }

  template< class Traits, class T, class A >
  inline OutStreamInterface< Traits > &
    operator<< ( OutStreamInterface< Traits > &out,
                 const std::vector< T, A > & value )
  {
    const size_t size = value.size();
    out << size;
    for( size_t i = 0; i < size; ++i )
      out << value[ i ];
    return out;
  }

  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 double &value )
  {
    in.readDouble( value );
    return in;
  }
  
  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 float &value )
  {
    in.readFloat( value );
    return in;
  }
  
  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 int &value )
  {
    in.readInt( value );
    return in;
  }
  
  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 char &value )
  {
    in.readChar( value );
    return in;
  }

  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 bool &value )
  {
    in.readBool( value );
    return in;
  }
  
  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 std :: string &s )
  {
    in.readString( s );
    return in;
  }
  
  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 unsigned int &value )
  {
    in.readUnsignedInt( value );
    return in;
  }

  template< class Traits >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 unsigned long &value )
  {
    in.readUnsignedLong( value );
    return in;
  }

  template< class Traits, class T, class A >
  inline InStreamInterface< Traits > &
    operator>> ( InStreamInterface< Traits > &in,
                 std::vector< T, A > & value )
  {
    size_t size = 0;
    in >> size;    
    value.resize( size );
    for( size_t i = 0; i < size; ++i )
      in >> value[ i ];
    return in;
  }

}

#endif // #ifndef DUNE_FEM_STREAMS_INLINE_HH
