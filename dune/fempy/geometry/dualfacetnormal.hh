#ifndef DUNE_FEM_MISC_DUALFACENORMAL_HH
#define DUNE_FEM_MISC_DUALFACENORMAL_HH

#include <dune/fem/quadrature.hh>

namespace Dune
{
  namespace Fem
  {
    namespace detail
    {
      template <class Geometry, class Point>
      static int
      point2Index( const Point& x )
      {
        // compute which point we are looking at
        // DualFacetQuadrature ;//
        return 0; //pt.index();
      }
      template <class Quadrature>
      static typename Quadrature::CoordinateType
      point2Index( const QuadraturePointWrapper< Quadrature >& pt)
      {
        return pt.index();
      }
    }


    template <class Entity, class Geometry, class Point>
    static
    typename Entity::Geometry::GlobalCoordinate
    dualFaceNormal( const Entity& entity,
                    const Geometry& geometry,
                    const Point& x,
                    const int face )
    {
      if constexpr ( Entity::dimension == 2 )
      {
        const auto& refElem = Dune::referenceElement( entity );
        const auto elemCenter = geometry.global( refElem.position(0,0) );
        const auto faceCenter = geometry.global( refElem.position(face,1) );

        auto normal = faceCenter - elemCenter;
        auto swap = normal[ 0 ];
        normal[ 0 ] = -normal[ 1 ];
        normal[ 1 ] = swap;

        // this need to be scaled since in GalerkinScheme the
        // geometry.integrationElement is multiplied which it should not in this
        // case.
        normal *= 1./geometry.volume();
        return normal;
      }
      else // dim == 3
      {
        // TODO: implement 3d version
        std::Abort();
      }

    }

    template <class Entity, class Geometry, class Point>
    static
    typename Entity::Geometry::GlobalCoordinate
    dualFaceNormal( const Entity& entity,
                    const Geometry& geometry,
                    const Point& x )
    {
      return dualFaceNormal( entity, geometry, x, point2Index( x ) );
    }

  } // end namespace Fem

} // end namespace Dune

#endif // DUNE_FEM_MISC_DUALFACENORMAL_HH
