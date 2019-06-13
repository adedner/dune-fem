#ifndef DUNE_FEMPY_QUADRATURE_CACHINGPOINT_HH
#define DUNE_FEMPY_QUADRATURE_CACHINGPOINT_HH

#include <type_traits>

#include <dune/common/ftraits.hh>
#include <dune/common/fvector.hh>

#include <dune/fem/quadrature/cachingpointlist.hh>

namespace Dune
{

  namespace FemPy
  {

    // CachingPoint
    // ------------

    template< class Coordinate, int codim = 0 >
    struct CachingPoint;

    template< class Coordinate >
    struct CachingPoint< Coordinate, 0 >
      : public Fem::CachingInterface
    {
      typedef CachingPoint< Coordinate, 0 > This;

    public:
      typedef Coordinate CoordinateType;
      typedef typename FieldTraits< Coordinate >::real_type RealType;
      typedef Coordinate LocalCoordinateType;

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      CachingPoint ( const Quadrature &quadrature, std::size_t idx )
        : id_( quadrature.id() ),
          idx_( quadrature.cachingPoint( idx ) ),
          position_( quadrature.point( idx ) )
      {}

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      explicit CachingPoint ( const Fem::QuadraturePointWrapper< Quadrature > &x )
        : CachingPoint( x.quadrature(), x.index() )
      {}

      explicit operator Fem::QuadraturePointWrapper< This > () const noexcept { return Fem::QuadraturePointWrapper< This >( *this, 0u ); }

      std::size_t id () const { return id_; }

      const CoordinateType &point ( std::size_t qp ) const { return position_; }
      const LocalCoordinateType &localPoint ( std::size_t qp ) const { return position_; }

      std::size_t cachingPoint ( std::size_t qp ) const { assert( qp == 0u ); return idx_; }

    private:
      std::size_t id_, idx_;
      const CoordinateType &position_;
    };

    template< class Coordinate >
    struct CachingPoint< Coordinate, 1 >
      : public Fem::CachingInterface
    {
      typedef CachingPoint< Coordinate, 1 > This;

    public:
      typedef Coordinate CoordinateType;
      typedef typename FieldTraits< Coordinate >::real_type RealType;
      typedef FieldVector< typename FieldTraits< Coordinate >::field_type, Coordinate::dimension-1 > LocalCoordinateType;

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      CachingPoint ( const Quadrature &quadrature, std::size_t idx )
        : id_( quadrature.id() ),
          idx_( quadrature.cachingPoint( idx ) ),
          position_( quadrature.point( idx ) ),
          localPosition_( quadrature.localPoint( idx ) )
      {}

      template< class Quadrature, std::enable_if_t< std::is_convertible< Quadrature, Fem::CachingInterface >::value, int > = 0 >
      explicit CachingPoint ( const Fem::QuadraturePointWrapper< Quadrature > &x )
        : CachingPoint( x.quadrature(), x.index() )
      {}

      explicit operator Fem::QuadraturePointWrapper< This > () const noexcept { return Fem::QuadraturePointWrapper< This >( *this, 0u ); }

      std::size_t id () const { return id_; }

      const CoordinateType &point ( std::size_t qp ) const { return position_; }
      const LocalCoordinateType &localPoint ( std::size_t qp ) const { return localPosition_; }

      std::size_t cachingPoint ( std::size_t qp ) const { assert( qp == 0u ); return idx_; }

    private:
      std::size_t id_, idx_;
      const CoordinateType &position_;
      const LocalCoordinateType &localPosition_;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_QUADRATURE_CACHINGPOINT_HH
