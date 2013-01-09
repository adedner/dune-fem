#ifndef DUNE_FEM_BASEFUNCTIONSET_FUNCTOR_HH
#define DUNE_FEM_BASEFUNCTIONSET_FUNCTOR_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/misc/functor.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class T >
    inline void axpy ( const T &a, const T &x, T &y );

    template< class K, int SIZE >
    inline void axpy ( const typename FieldTraits< K >::field_type &a,
                       const FieldVector< K, SIZE > &x,
                       FieldVector< K, SIZE > &y );

    template< class K, int ROWS, int COLS >
    inline void axpy ( const typename FieldTraits< K >::field_type &a,
                       const FieldMatrix< K, ROWS, COLS > &x,
                       FieldMatrix< K, ROWS, COLS > &y );



    // axpy
    // ----

    template< class T >
    inline void axpy ( const T &a, const T &x, T &y )
    {
      y += a*x;
    }

    template< class K, int SIZE >
    inline void axpy ( const typename FieldTraits< K >::field_type &a,
                       const FieldVector< K, SIZE > &x,
                       FieldVector< K, SIZE > &y )
    {
      for( int i = 0; i < SIZE; ++i )
        axpy( a, x[ i ], y[ i ] );
    }

    template< class K, int ROWS, int COLS >
    inline void axpy ( const typename FieldTraits< K >::field_type &a,
                       const FieldMatrix< K, ROWS, COLS > &x,
                       FieldMatrix< K, ROWS, COLS > &y )
    {
      y.axpy( a, x );
    }



    // scalarProduct
    // -------------

    template< class T >
    inline typename T::field_type scalarProduct ( const T &a, const T &b )
    {
      return a * b;
    }

    template< class K, int ROWS, int COLS >
    inline K scalarProduct ( const FieldMatrix< K, ROWS, COLS > &a, const FieldMatrix< K, ROWS, COLS > &b )
    {
      K s( 0 );
      for( int r = 0; r < ROWS; ++r )
        s += a[ r ] * b[ r ];
      return s;
    }



    // AxpyFunctor
    // -----------

    template< class Vector, class Value >
    struct AxpyFunctor
    {
      AxpyFunctor ( const Vector &vector, Value &value )
      : vector_( vector ),
        value_( value )
      {}

      void operator() ( const std::size_t i, const Value &v )
      {
        axpy( vector_[ i ], v, value_ );
      }

    private:
      const Vector &vector_;
      Value &value_;
    };



    // FunctionalAxpyFunctor
    // ---------------------

    template< class Value, class Vector >
    struct FunctionalAxpyFunctor
    {
      FunctionalAxpyFunctor ( const Value &value, Vector &vector )
      : value_( value ),
        vector_( vector )
      {}

      void operator() ( const std::size_t i, const Value &v )
      {
        vector_[ i ] += scalarProduct( v, value_ );
      }

    private:
      const Value &value_;
      Vector &vector_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASEFUNCTIONSET_FUNCTOR_HH
