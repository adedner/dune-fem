#ifndef DUNE_FEM_SPACE_FOURIER_BASISFUNCTIONS_HH
#define DUNE_FEM_SPACE_FOURIER_BASISFUNCTIONS_HH

#include <cassert>
#include <cstddef>
#include <limits>

#include <dune/common/array.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/power.hh>
#include <dune/common/static_assert.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    // NumFourierBasisFunctions
    // ------------------------

    template< int dimension, int order >
    struct NumFourierBasisFunctions
    {
      static const int v = StaticPower< (2*order+1), dimension >::power;

    };



    // FourierBasisFunctions
    // ---------------------

    template< class FunctionSpace, int order >
    struct FourierBasisFunctions
    {
      dune_static_assert( (FunctionSpace::dimRange == 1),
                          "FunctionSpace must be scalar (i.e., dimRange = 1)." );

      typedef FunctionSpace FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef std::size_t SizeType;

      //! \brief return number of basis functions
      SizeType size () const;

      /**
       * \brief evalute each basis function
       *
       *  \param[in]  x        global coordinate
       *  \param[in]  functor  functor call for evaluating each basis function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
       *  \code
       *  struct Functor
       *  {
       *    template< class Value >
       *    void operator() ( const int basisFunction, const Value &value );
       *  };
       *  \endcode
       */
      template< class Functor >
      static void evaluateEach ( const DomainType &x, Functor functor );

      /**
       * \brief evalute jacobian of each basis function
       *
       *  \param[in]  x        global coordinate
       *  \param[in]  functor  functor call for evaluating the jacobian of each basis function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
       *  \code
       *  struct Functor
       *  {
       *    template< class Jacobian >
       *    void operator() ( const int basisFunction, const Jacobian &jacobian );
       *  };
       *  \endcode
       */
      template< class Functor >
      static void jacobianEach ( const DomainType &x, Functor functor );

      /**
       * \brief evalute hessian of each basis function
       *
       *  \param[in]  x        global coordinate
       *  \param[in]  functor  functor call for evaluating the hessian of each basis function
       *
       *  The functor has to be a copyable object satisfying the following
       *  interface:
       *  \code
       *  struct Functor
       *  {
       *    template< class Hessian >
       *    void operator() ( const int basisFunction, const Hessian &hessian );
       *  };
       *  \endcode
       */
      template< class Functor >
      static void hessianEach ( const DomainType &x, Functor functor );
    };



    // Template specialization for dimDomain = dimRange = 1
    // ----------------------------------------------------

    template< class DomainFieldType, class RangeFieldType, int order >
    class FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, order >
    {
      typedef FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, order > ThisType;

    public:
      typedef FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 > FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef std::size_t SizeType;

      static SizeType size () { return NumFourierBasisFunctions< 1, order >::v; }

      template< class Functor >
      static void evaluateEach ( const DomainType &x, Functor functor )
      {
        functor( 0, .5 );
        for( SizeType n = 1; n <= order; ++n )
        {
          const int basisFunction = 2*n-1;
          functor( basisFunction, cos( n*x ) );
          functor( basisFunction+1, sin( n*x ) );
        }
      }

      template< class Functor >
      static void jacobianEach ( const DomainType &x, Functor functor )
      {
        functor( 0, 0. );
        for( SizeType n = 1; n <= order; ++n )
        {
          const int basisFunction = 2*n-1;
          functor( basisFunction, -n*sin( n*x ) );
          functor( basisFunction+1, n*cos( n*x ) );
        }
      }

      template< class Functor >
      static void hessianEach ( const DomainType &x, Functor functor )
      {
        functor( 0, 0. );
        for( SizeType n = 1; n <= order; ++n )
        {
          const int basisFunction = 2*n-1;
          functor( basisFunction, -n*n*cos( n*x ) );
          functor( basisFunction+1, -n*n*sin( n*x ) );
        }
      }
    };



    // Template specialization for dimDomain > 1, dimRange = 1
    // -------------------------------------------------------

    template< class DomainFieldType, class RangeFieldType, int dimDomain, int order >
    class FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, order >
    {
      typedef FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, order > ThisType;

    public:
      typedef FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 > FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef std::size_t SizeType;

    private:
      // number of Fourier basis function for dimDomain = 1
      static const int buffer_size = NumFourierBasisFunctions< 1, order >::v;

      // tags used for building cache
      struct Evaluate {};  //< evaluate basis functions
      struct Jacobian {};  //< evaluate basis functions and jacobians
      struct Hessian {};   //< evaluate basis functions, jacobians, and hessians

      struct Assign;

    protected:
      // multi index used for tensor product basis functions
      typedef Dune::FieldVector< int, dimDomain > MultiIndexType;

      // iterator type and methods for accessing multi indices
      struct MultiIndexIterator;
      typedef MultiIndexIterator IteratorType;
      static IteratorType begin () { return IteratorType::begin(); }
      static IteratorType end () { return IteratorType::end(); }

    public:
      static SizeType size () { return NumFourierBasisFunctions< dimDomain, order >::v; }

      template< class Functor >
      void evaluateEach ( const DomainType &x, Functor functor ) const
      {
        prepare( Evaluate(), x );
        SizeType index( 0 );
        const IteratorType end = ThisType::end();
        for( IteratorType it = ThisType::begin(); it != end; ++it, ++index )
        {
          RangeType value;
          evaluate( *it, value );
          assert( index == IteratorType::index( *it ) );
          functor( index, value );
        }
        assert( index == size() );
      }

      template< class Functor >
      void jacobianEach ( const DomainType &x, Functor functor ) const
      {
        prepare( Jacobian(), x );
        SizeType index( 0 );
        const IteratorType end = ThisType::end();
        for( IteratorType it = ThisType::begin(); it != end; ++it, ++index )
        {
          JacobianRangeType jacobian;
          evaluate( *it, jacobian );
          assert( index == IteratorType::index( *it ) );
          functor( index, jacobian );
        }
        assert( index == size() );
      }

      template< class Functor >
      void hessianEach ( const DomainType &x, Functor functor ) const
      {
        prepare( Hessian(), x );
        SizeType index( 0 );
        const IteratorType end = ThisType::end();
        for( IteratorType it = ThisType::begin(); it != end; ++it, ++index )
        {
          HessianRangeType hessian;
          evaluate( *it, hessian );
          assert( index == IteratorType::index( *it ) );
          functor( index, hessian );
        }
        assert( index == size() );
      }

    protected:
      // evaluate tensor product basis function
      void evaluate ( const MultiIndexType &multiIndex, RangeType &value ) const
      {
        value = RangeType( RangeFieldType( 1 ) );
        for( SizeType i = 0; i < dimDomain; ++i )
          value *= buffer_[ i ][ multiIndex[ i ] ];
      }

      // evaluate jacobian of tensor product basis function
      void evaluate ( const MultiIndexType &multiIndex, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( 1 );
        for( SizeType k = 0; k < dimDomain; ++k )
        {
          const RangeFieldType phi = buffer_[ k ][ multiIndex[ k ] ];
          const RangeFieldType dphi = buffer_[ k ][ buffer_size + multiIndex[ k ] ];
          for( int i = 0; i < dimDomain; ++i )
            jacobian[ 0 ][ i ] *= ( k == i ) ? dphi : phi;
        }
      }

      // evaluate hessian of tensor product basis function
      void evaluate ( const MultiIndexType &multiIndex, HessianRangeType &hessian ) const
      {
        hessian = HessianRangeType( 1 );
        for( int k = 0; k < dimDomain; ++k )
        {
          const RangeFieldType phi = buffer_[ k ][ multiIndex[ k ] ];
          const RangeFieldType dphi = buffer_[ k ][ buffer_size + multiIndex[ k ] ];
          for( int i = 0; i < dimDomain; ++i )
          {
            hessian[ i ][ i ] *= ( k == i ) ? buffer_[ i ][ 2*buffer_size + multiIndex[ i ] ] : phi;
            for( int j = i+1; j < dimDomain; ++j )
            {
              RangeFieldType tmp = ( k == i || k == j ) ? dphi : phi;
              hessian[ i ][ j ] *= tmp;
              hessian[ j ][ i ] *= tmp;
            }
          }
        }
      }

      // methods for building cache
      void prepare ( const Evaluate &, const DomainType &x ) const;
      void prepare ( const Jacobian &, const DomainType &x ) const;
      void prepare ( const Hessian &, const DomainType &x ) const;

    private:
      // cache for evaluation of basis functions
      mutable array< array< RangeFieldType, 3*buffer_size>, dimDomain > buffer_;
    };



    // Implementation of FourierBasisFunctions::Assign
    // -----------------------------------------------

    template< class DomainFieldType, class RangeFieldType, int dimDomain, int order >
    struct FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, order >::Assign
    {
      explicit Assign ( RangeFieldType *buffer ) : buffer_( buffer ) {}

      void operator() ( const std::size_t i, const RangeFieldType &value )
      {
        buffer_[ i ] = value;
      }

      template< class T >
      void operator() ( const std::size_t i, const FieldVector< T, 1 > &value )
      {
        (*this)( i, value[ 0 ] );
      }

      template< class T >
      void operator() ( const std::size_t i, const FieldMatrix< T, 1, 1 > &value )
      {
        (*this)( i, value[ 0 ][ 0 ] );
      }

    private:
      RangeFieldType *buffer_;
    };



    // Implementation of FourierBasisFunctions::MultiIndexIterator
    // -----------------------------------------------------------

    template< class DomainFieldType, class RangeFieldType, int dimDomain, int order >
    struct FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, order >::MultiIndexIterator
    {
      typedef MultiIndexIterator ThisType;

    protected:
      typedef int IndexType;

      explicit MultiIndexIterator ( IndexType n ) : multiIndex_( n ) {}

      static IndexType invalidIndex () { return std::numeric_limits< IndexType >::max(); }

      static const int N = NumFourierBasisFunctions< 1, order >::v;

    public:
      static ThisType begin () { return ThisType( 0 ); }
      static ThisType end () { return ThisType( invalidIndex() ); }

      ThisType operator++ ()
      {
        // try to increment and leave eventually
        for( int i = 0; i < dimDomain; ++i )
        {
          const int j = dimDomain-i-1;
          if( ++multiIndex_[ j ] <= N )
            return *this;
          multiIndex_[ j ] = 0;
        }
        
        // otherwise, reset this iterator to end iterator
        *this = end(); 
        return *this;
      }

      bool operator== ( const ThisType &other ) const { return ( multiIndex_ == other.multiIndex_ ); }

      bool operator!= ( const ThisType &other ) const { return !( *this == other ); }

      const MultiIndexType &operator* () const { return multiIndex_; }

      const MultiIndexType *operator-> () const { return &multiIndex_; }

      SizeType index () const { return index( multiIndex_ ); }

      static SizeType index ( const MultiIndexType &multiIndex )
      {
        SizeType index = 0, factor = 1;
        for( SizeType i = dimDomain-1; i > 0; --i )
        {
          index += multiIndex[ i ]*factor;
          factor *= N+1;
        }
        index += multiIndex[ 0 ]*factor;
        assert( index < size() );
        return index;
      }

    private:
      MultiIndexType multiIndex_;
    };



    // Implementation of FourierBasisFunctions
    // ---------------------------------------

    template< class DomainFieldType, class RangeFieldType, int dimDomain, int order >
    void FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, order >
      ::prepare ( const Evaluate &, const DomainType &x ) const
    {
      typedef FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, order > BasisFunctionsImp;
      for( SizeType i = 0; i < dimDomain; ++i )
      {
        RangeFieldType *it = &buffer_[ i ][ 0 ];
        Dune::FieldVector< DomainFieldType, 1 > y( x[ i ] );
        BasisFunctionsImp::evaluateEach( y, Assign( it ) );
      }
    };


    template< class DomainFieldType, class RangeFieldType, int dimDomain, int order >
    void FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, order >
      ::prepare ( const Jacobian &, const DomainType &x ) const
    {
      typedef FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, order > BasisFunctionsImp;
      for( SizeType i = 0; i < dimDomain; ++i )
      {
        RangeFieldType *it = &buffer_[ i ][ 0 ];
        Dune::FieldVector< DomainFieldType, 1 > y( x[ i ] );
        BasisFunctionsImp::evaluateEach( y, Assign( it ) );
        BasisFunctionsImp::jacobianEach( y, Assign( it+buffer_size ) );
      }
    };


    template< class DomainFieldType, class RangeFieldType, int dimDomain, int order >
    void FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, 1 >, order >
      ::prepare ( const Hessian &, const DomainType &x ) const
    {
      typedef FourierBasisFunctions< FunctionSpace< DomainFieldType, RangeFieldType, 1, 1 >, order > BasisFunctionsImp;
      for( SizeType i = 0; i < dimDomain; ++i )
      {
        RangeFieldType *it = &buffer_[ i ][ 0 ];
        Dune::FieldVector< DomainFieldType, 1 > y( x[ i ] );
        BasisFunctionsImp::evaluateEach( y, Assign( it ) );
        BasisFunctionsImp::jacobianEach( y, Assign( it+buffer_size) );
        BasisFunctionsImp::hessianEach( y, Assign( it+2*buffer_size ) );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_BASISFUNCTIONS_HH
