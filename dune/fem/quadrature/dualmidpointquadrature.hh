#ifndef DUNE_FEM_DUALMIDPOINTQUADRATURE_HH
#define DUNE_FEM_DUALMIDPOINTQUADRATURE_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/lagrange/lagrangepoints.hh>

namespace Dune {

namespace Fem {

/**Define a lumping quadrature for all geometries. Note, however, that
 * this may not make sense for anything else than simplices or maybe
 * hexagonal grids. For simplicial meshes the quadrature formula is
 * exact on linear polynomials and hence the quadrature error is
 * quadratic in the mesh-size. A mass-matrix assembled with the
 * caching quadrature will be diagonal in the context of Lagrange
 * spaces. Generally, it is a bad idea to use mass-lumping for
 * anything else than linear (or maybe bilinear) finite elements.
 *
 * There are probably much more efficient ways to perform
 * mass-lumping. This "quadrature" approach is convenient, however, as
 * it can simply be plugged into existing code by replacing
 * the quadrature rules.
 */


template<class FieldImp, Dune::GeometryType::Id geometryId>
class DualFacetMidpointQuadrature
  : public QuadratureImp<FieldImp, Dune::GeometryType(geometryId).dim()>
{
public:
  typedef FieldImp FieldType;
  static constexpr auto dimension = Dune::GeometryType(geometryId).dim();

private:
  typedef DualFacetMidpointQuadrature<FieldType, geometryId> ThisType;
  typedef QuadratureImp<FieldType, dimension> BaseType;

  static const int maxOrder_ = 6;
  static const int order_ = 1;

public:
  typedef typename BaseType::CoordinateType CoordinateType;

  /** \brief constructor filling the list of points and weights.
   *
   *  \param[in]  gt     geometry type for which a quadrature is desired
   *  \param[in]  order  order (between 1 and 6)
   *  \param[in]  id     unique identifier, ignored
   */
  DualFacetMidpointQuadrature(const GeometryType& gt, int order, int id)
    : BaseType(id)
  {
    // create Lagrange points to define quadrature
    //typedef LagrangePointListImplementation< FieldType, Dune::GeometryType(geometryId).id(), dimension, maxOrder_ > LP;
    //LP pts( gt, order_, id );

    //const auto &refElement = Dune::ReferenceElements< FieldType, dimension >::general( gt );
    //const std::size_t numCorners = refElement.size( dimension );
    //const std::size_t numCorners = pts.nop();
    if( gt.isTriangle() )
    {
      CoordinateType p0( {{ (0.5+1./3.)*0.5, (1./6.) }} );
      this->addQuadraturePoint( p0, 1./3. );
      CoordinateType p1( {{ (1./6.), (0.5+1./3.)*0.5 }} );
      this->addQuadraturePoint( p1, 1./3. );
      CoordinateType p2( {{ (0.5+1./3.)*0.5, (0.5+1./3.)*0.5 }} );
      this->addQuadraturePoint( p2, 1./3. );
    }
    else if ( gt.isQuadrilateral() )
    {
      CoordinateType p0( {{ 0.25, 0.5 }} );
      this->addQuadraturePoint( p0, 0.25 );
      CoordinateType p1( {{ 0.75, 0.5 }} );
      this->addQuadraturePoint( p1, 0.25 );

      CoordinateType p2( {{ 0.5, 0.25 }} );
      this->addQuadraturePoint( p2, 0.25 );
      CoordinateType p3( {{ 0.5, 0.75 }} );
      this->addQuadraturePoint( p2, 0.25 );
    }
    else
    {
      DUNE_THROW(NotImplemented,"not implemented yet!");
    }
  }

  /** \copydoc QuadratureImp::geometry
   */
  virtual GeometryType geometryType() const { return Dune::GeometryType(geometryId); }
  /** \copydoc QuadratureImp::order
   */
  virtual int order () const { return order_; }

  //! maximal order of available quadratures
  static std::size_t maxOrder () { return maxOrder_; }
};

template<class FieldType, int dimension>
struct DefaultDualFacetMidpointQuadratureTraits
{
  typedef QuadratureImp<FieldType, dimension> IntegrationPointListType;

  static constexpr Dune::GeometryType::Id simplexId = Dune::GeometryTypes::simplex(dimension);
  static constexpr Dune::GeometryType::Id cubeId    = Dune::GeometryTypes::cube(dimension);
  static constexpr Dune::GeometryType::Id prismId   = Dune::GeometryTypes::prism ;
  static constexpr Dune::GeometryType::Id pyramidId = Dune::GeometryTypes::pyramid;

  typedef DualFacetMidpointQuadrature<FieldType, simplexId>  SimplexQuadratureType;
  typedef DualFacetMidpointQuadrature<FieldType, cubeId   >  CubeQuadratureType;
  typedef DualFacetMidpointQuadrature<FieldType, prismId  >  PrismQuadratureType;
  typedef DualFacetMidpointQuadrature<FieldType, pyramidId>  PyramidQuadratureType;

  typedef SimplexQuadratureType LineQuadratureType;
  typedef SimplexQuadratureType PointQuadratureType;

  typedef int QuadratureKeyType;
};

// DualFacetMidpointQuadrature uses CachingQuadrature with a different traits class for creating the quadratures.
template<class GridPartImp, int codim>
using CachingDualFacetMidpointQuadrature = CachingQuadrature< GridPartImp, codim, DefaultDualFacetMidpointQuadratureTraits >;

} // Fem

} // Dune


#endif // header guard
