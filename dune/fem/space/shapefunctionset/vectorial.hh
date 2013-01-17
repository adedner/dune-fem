#ifndef DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH

// C++ includes
#include <algorithm>
#include <cstddef>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>


namespace Dune
{

  namespace Fem 
  {

    // MakeVectorialTraits
    // -------------------

    template< class Scalar, class Vectorial >
    struct MakeVectorialTraits;

    template< class K, int dimR >
    struct MakeVectorialTraits< FieldVector< K, 1 >, FieldVector< K, dimR > >
    {
      typedef FieldVector< K, 1 > ScalarType;
      typedef FieldVector< K, dimR > VectorialType;

      typedef typename FieldTraits< VectorialType >::field_type field_type;
      typedef typename VectorialType::size_type ComponentType;
      typedef typename VectorialType::size_type size_type;

      static const size_type factor = dimR;

      static ComponentType begin () { return ComponentType( 0 ); }
      static ComponentType end () { return ComponentType( factor ); }

      static VectorialType zeroVectorial () { return VectorialType( K( field_type( 0 ) ) ); }

      static const K &access ( const ScalarType &x ) { return x[ 0 ]; }
      static K &access ( ScalarType &x ) { return x[ 0 ]; }

      static const K &access ( const VectorialType &x, const ComponentType &i ) { return x[ i ]; }
      static K &access ( VectorialType &x, const ComponentType &i ) { return x[ i ]; }

      static size_type index ( const ComponentType &i ) { return i; }
    };

    template< class K, int dimR, int dimD >
    struct MakeVectorialTraits< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > >
    {
      typedef FieldMatrix< K, 1, dimD > ScalarType;
      typedef FieldMatrix< K, dimR, dimD > VectorialType;

      typedef typename FieldTraits< VectorialType >::field_type field_type;
      typedef typename VectorialType::size_type ComponentType;
      typedef typename VectorialType::size_type size_type;

      static const size_type factor = dimR;

      static ComponentType begin () { return ComponentType( 0 ); }
      static ComponentType end () { return ComponentType( factor ); }

      static VectorialType zeroVectorial () { return VectorialType( K( field_type( 0 ) ) ); }

      static const FieldVector< K, dimD > &access ( const ScalarType &x ) { return x[ 0 ]; }
      static FieldVector< K, dimD > &access ( ScalarType &x ) { return x[ 0 ]; }

      static const FieldVector< K, dimD > &access ( const VectorialType &x, const ComponentType &i ) { return x[ i ]; }
      static FieldVector< K, dimD > &access ( VectorialType &x, const ComponentType &i ) { return x[ i ]; }

      static size_type index ( const ComponentType &i ) { return i; }
    };



    // MakeVectorialExpression
    // -----------------------

    template< class Scalar, class Vectorial >
    class BasicMakeVectorialExpression
    {
      typedef BasicMakeVectorialExpression< Scalar, Vectorial > ThisType;

      typedef MakeVectorialTraits< Scalar, Vectorial > Traits;

    public:
      typedef typename Traits::ScalarType ScalarType;
      typedef typename Traits::VectorialType VectorialType;

      typedef typename Traits::field_type field_type;
      typedef typename Traits::ComponentType ComponentType;
      typedef typename Traits::size_type size_type;

      BasicMakeVectorialExpression ( const ComponentType &component, const ScalarType &scalar )
      : component_( component ),
        scalar_( scalar )
      {}

      operator VectorialType () const
      {
        VectorialType vectorial = Traits::zeroVectorial();
        Traits::access( vectorial, component() ) = Traits::access( scalar() );
        return vectorial;
      }

      const ThisType &operator*= ( const field_type &s )
      {
        scalar() *= s;
        return *this;
      }

      const ThisType &operator/= ( const field_type &s )
      {
        scalar() /= s;
        return *this;
      }

      const ComponentType &component () const { return component_; }

      const ScalarType &scalar () const { return scalar_; }
      ScalarType &scalar () { return scalar_; }

    protected:
      ComponentType component_;
      ScalarType scalar_;
    };



    // MakeVectorialExpression
    // -----------------------

    template< class Scalar, class Vectorial >
    class MakeVectorialExpression
    : public BasicMakeVectorialExpression< Scalar, Vectorial >
    {
      typedef MakeVectorialExpression< Scalar, Vectorial > ThisType;
      typedef BasicMakeVectorialExpression< Scalar, Vectorial > BaseType;
      
    public:
      typedef typename BaseType::ComponentType ComponentType;
      typedef typename BaseType::ScalarType ScalarType;

      MakeVectorialExpression ( const ComponentType &component, const ScalarType &scalar )
      : BaseType( component, scalar )
      {}
    };

    template< class K, int dimR >
    class MakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > >
    : public BasicMakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > >
    {
      typedef MakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > > ThisType;
      typedef BasicMakeVectorialExpression< FieldVector< K, 1 >, FieldVector< K, dimR > > BaseType;
      
    public:
      typedef typename BaseType::ScalarType ScalarType;
      typedef typename BaseType::VectorialType VectorialType;

      typedef typename BaseType::field_type field_type;
      typedef typename BaseType::ComponentType ComponentType;
      typedef typename BaseType::size_type size_type;

      using BaseType::component;
      using BaseType::scalar;

      MakeVectorialExpression ( const ComponentType &component, const ScalarType &scalar )
      : BaseType( component, scalar )
      {}

      field_type operator* ( const ThisType &other ) const
      {
        return (component() == other.component() ? scalar() * other.scalar() : field_type( 0 ));
      }

      field_type operator* ( const VectorialType &other ) const
      {
        return (scalar()[ 0 ] * other[ component() ]);
      }

      field_type one_norm () const { return scalar().one_norm(); }
      field_type two_norm () const { return scalar().two_norm(); }
      field_type two_norm2 () const { return scalar().two_norm2(); }
      field_type infinity_norm () const { return scalar().infinity_norm(); }

      size_type size () const { return dimR; }

      friend field_type operator* ( const VectorialType &a, ThisType &b ) { return b*a; }
    };

    template< class K, int dimR, int dimD >
    class MakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > >
    : public BasicMakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > >
    {
      typedef MakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > > ThisType;
      typedef BasicMakeVectorialExpression< FieldMatrix< K, 1, dimD >, FieldMatrix< K, dimR, dimD > > BaseType;
      
    public:
      typedef typename BaseType::ScalarType ScalarType;
      typedef typename BaseType::VectorialType VectorialType;

      typedef typename BaseType::field_type field_type;
      typedef typename BaseType::ComponentType ComponentType;
      typedef typename BaseType::size_type size_type;

      using BaseType::component;
      using BaseType::scalar;

      MakeVectorialExpression ( const ComponentType &component, const ScalarType &scalar )
      : BaseType( component, scalar )
      {}

      template< class X, class Y >
      void mv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] = field_type( 0 );
        for( size_type j= 0; j < cols(); ++j )
          y[ component() ] += scalar()[ component() ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void mtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] = scalar()[ i ][ component() ] * x[ component() ];
      }

      template< class X, class Y >
      void umv ( const X &x, Y &y ) const
      {
        for( size_type j= 0; j < cols(); ++j )
          y[ component() ] += scalar()[ component() ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void umtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] += scalar()[ i ][ component() ] * x[ component() ];
      }

      template< class X, class Y >
      void mmv ( const X &x, Y &y ) const
      {
        for( size_type j= 0; j < cols(); ++j )
          y[ component() ] -= scalar()[ component() ][ j ] * x[ j ];
      }

      template< class X, class Y >
      void mmtv ( const X &x, Y &y ) const
      {
        for( size_type i = 0; i < rows(); ++i )
          y[ i ] -= scalar()[ i ][ component() ] * x[ component() ];
      }

      field_type frobenius_norm () const { return scalar().frobenius_norm(); }
      field_type frobenius_norm2 () const { return scalar().frobenius_norm2(); }
      field_type infinity_norm () const { return scalar().infinity_norm(); }

      field_type determinant () const { return (dimR == 1 ? scalar().determinant() : field_type( 0 )); }

      size_type N () const { return rows(); }
      size_type M () const { return cols(); }

      size_type rows () const { return dimR; }
      size_type cols () const { return scalar().cols(); }
    };



    // Auxilliary Functions for MakeVectorialExpression
    // ------------------------------------------------

    template< class Scalar, class Vectorial >
    inline bool
    operator== ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                 const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return ((a.component() == b.component()) && (a.scalar() == b.scalar()));
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator!= ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                 const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return ((a.component() != b.component()) || (a.scalar() != b.scalar()));
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator== ( const Vectorial &a,
                 const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return (a == static_cast< Vectorial >( b ));
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator!= ( const Vectorial &a,
                 const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return (a != static_cast< Vectorial >( b ));
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator== ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                 const Vectorial &b )
    {
      return (static_cast< Vectorial >( a ) == b);
    }

    template< class Scalar, class Vectorial >
    inline bool
    operator!= ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                 const Vectorial &b )
    {
      return (static_cast< Vectorial >( a ) != b);
    }

    template< class Scalar, class Vectorial >
    inline void
    axpy ( const typename MakeVectorialTraits< Scalar, Vectorial >::field_type &a,
           const MakeVectorialExpression< Scalar, Vectorial > &x,
           typename MakeVectorialTraits< Scalar, Vectorial >::VectorialType &y )
    {
      typedef MakeVectorialTraits< Scalar, Vectorial > Traits;
      axpy( a, Traits::access( x.scalar() ), Traits::access( y, x.component() ) );
    }

    template< class GeometryJacobianInverseTransposed, class K, int ROWS >
    void jacobianTransformation ( const GeometryJacobianInverseTransposed &gjit,
                                  const MakeVectorialExpression< FieldMatrix< K, 1, GeometryJacobianInverseTransposed::rows >, FieldMatrix< K, ROWS, GeometryJacobianInverseTransposed::cols > > &a,
                                  FieldMatrix< K, ROWS, GeometryJacobianInverseTransposed::rows > &b )
    {
      typedef MakeVectorialTraits< FieldMatrix< K, 1, GeometryJacobianInverseTransposed::rows >, FieldMatrix< K, ROWS, GeometryJacobianInverseTransposed::cols > > Traits;
      b = Traits::zeroVectorial();
      gjit.mv( Traits::access( a.scalar() ), Traits::access( b, a.component() ) );
    }

    template< class GeometryJacobianInverseTransposed, class K, int SIZE >
    void hessianTransformation ( const GeometryJacobianInverseTransposed &gjit,
                                 const MakeVectorialExpression< FieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::cols, GeometryJacobianInverseTransposed::cols >, 1 >, FieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::cols, GeometryJacobianInverseTransposed::cols >, SIZE > > &a,
                                 FieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::rows, GeometryJacobianInverseTransposed::rows >, SIZE > &b )
    {
      typedef MakeVectorialTraits< FieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::cols, GeometryJacobianInverseTransposed::cols >, 1 >, FieldVector< FieldMatrix< K, GeometryJacobianInverseTransposed::cols, GeometryJacobianInverseTransposed::cols >, SIZE > > Traits;
      b = Traits::zeroVectorial();
      for( int k = 0; k < GeometryJacobianInverseTransposed::cols; ++k )
      {
        FieldVector< K, GeometryJacobianInverseTransposed::rows > c;
        gjit.mv( Traits::access( a.scalar() )[ k ], c );
        for( int j = 0; j < GeometryJacobianInverseTransposed::rows; ++j )
          Traits::access( b, a.component() )[ j ].axpy( gjit[ j ][ k ], c );
      }
    }

    template< class Scalar, class Vectorial >
    inline typename MakeVectorialTraits< Scalar, Vectorial >::field_type
    scalarProduct ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                    const Vectorial &b )
    {
      return scalarProduct( a.scalar()[ 0 ], b[ a.component() ] );
    }

    template< class Scalar, class Vectorial >
    inline typename MakeVectorialTraits< Scalar, Vectorial >::field_type
    scalarProduct ( const Vectorial &a,
                    const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      return scalarProduct( a[ b.component() ], b.scalar()[ 0 ] );
    }

    template< class Scalar, class Vectorial >
    inline typename MakeVectorialTraits< Scalar, Vectorial >::field_type
    scalarProduct ( const MakeVectorialExpression< Scalar, Vectorial > &a,
                    const MakeVectorialExpression< Scalar, Vectorial > &b )
    {
      typedef typename MakeVectorialTraits< Scalar, Vectorial >::field_type field_type;
      return (a.component() == b.component() ? scalarProduct( a.scalar(), b.scalar() ) : field_type( 0 ));
    }



    // ToNewRange
    // ----------

    template< class ScalarFunctionSpace, class RangeVector >
    struct ToNewRange;

    template< class DomainField, class RangeField, int dimD, int dimR >
    struct ToNewRange< FunctionSpace< DomainField, RangeField, dimD, 1 >, FieldVector< RangeField, dimR > >
    {
      typedef FunctionSpace< DomainField, RangeField, dimD, dimR > Type;
    };



    // VectorialShapeFunctionSet
    // -------------------------

    template< class ScalarShapeFunctionSet, class RangeVector >
    class VectorialShapeFunctionSet
    {
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector > ThisType;

    public:
      typedef ScalarShapeFunctionSet ScalarShapeFunctionSetType;

    protected:
      typedef typename ScalarShapeFunctionSetType::FunctionSpaceType ScalarFunctionSpaceType;

      static const std::size_t dimRangeFactor = MakeVectorialTraits< typename ScalarFunctionSpaceType::RangeType, RangeVector >::factor;

      template< class Functor, class Vectorial >
      struct VectorialFunctor;

    public:
      typedef typename ToNewRange< ScalarFunctionSpaceType, RangeVector >::Type FunctionSpaceType;

      explicit VectorialShapeFunctionSet ( const ScalarShapeFunctionSetType &scalarShapeFunctionSet = ScalarShapeFunctionSetType() )
      : scalarShapeFunctionSet_( scalarShapeFunctionSet )
      {}

      const ScalarShapeFunctionSetType &scalarShapeFunctionSet () const { return scalarShapeFunctionSet_; }

      // Shape Function Set Interface Methods
      std::size_t size () const { return dimRangeFactor * scalarShapeFunctionSet().size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const;
     
    protected:
      ScalarShapeFunctionSet scalarShapeFunctionSet_;
    };



    // VectorialShapeFunctionSet::VectorialFunctor
    // -------------------------------------------

    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Functor, class Vectorial >
    struct VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >::VectorialFunctor
    {
      explicit VectorialFunctor ( const Functor &functor )
      : functor_( functor )
      {}

      template< class Scalar >
      void operator() ( const std::size_t i, const Scalar &value )
      {
        typedef MakeVectorialTraits< Scalar, Vectorial > Traits;
        typedef MakeVectorialExpression< Scalar, Vectorial > Expression;
        const typename Traits::ComponentType end = Traits::end();
        for( typename Traits::ComponentType k = Traits::begin(); k != end; ++k )
          functor_( i*Traits::factor + Traits::index( k ), Expression( k, value ) );
      }

    private:
      Functor functor_;
    };



    // Implementation of VectorialShapeFunctionSet
    // -------------------------------------------

    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Point, class Functor >
    inline void VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >
      ::evaluateEach ( const Point &x, Functor functor ) const
    {
      typedef typename FunctionSpaceType::RangeType VectorialType;
      scalarShapeFunctionSet().evaluateEach( x, VectorialFunctor< Functor, VectorialType >( functor ) );
    }
    

    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Point, class Functor >
    inline void VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >
      ::jacobianEach ( const Point &x, Functor functor ) const
    {
      typedef typename FunctionSpaceType::JacobianRangeType VectorialType;
      scalarShapeFunctionSet().jacobianEach( x, VectorialFunctor< Functor, VectorialType >( functor ) );
    }


    template< class ScalarShapeFunctionSet, class RangeVector >
    template< class Point, class Functor >
    inline void VectorialShapeFunctionSet< ScalarShapeFunctionSet, RangeVector >
      ::hessianEach ( const Point &x, Functor functor ) const
    {
      typedef typename FunctionSpaceType::HessianRangeType VectorialType;
      scalarShapeFunctionSet().hessianEach( x, VectorialFunctor< Functor, VectorialType >( functor ) );
    }

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_VECTORIAL_HH
