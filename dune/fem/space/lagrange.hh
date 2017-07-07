#ifndef HAVE_DUNE_FEM_SPACE_LAGRANGE
#define HAVE_DUNE_FEM_SPACE_LAGRANGE

#include <dune/fem/space/lagrange/space.hh>

#if HAVE_DUNE_LOCALFUNCTIONS
// dune-localfunctions includes

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/localfiniteelement/space.hh>

namespace Dune
{
  namespace Fem
  {

    template< class FunctionSpace, class GridPart,
              template< class, unsigned int > class PointSet = EquidistantPointSet
            >
    class LagrangeFiniteElementMap
    {
      typedef LagrangeFiniteElementMap< FunctionSpace, GridPart > ThisType;

    public:
      typedef GridPart GridPartType;

      typedef unsigned int KeyType;

      typedef typename FunctionSpace::DomainFieldType DomainFieldType;
      typedef typename FunctionSpace::RangeFieldType RangeFieldType;

      static const int dimLocal = GridPart::dimension;

      typedef LagrangeLocalFiniteElement< PointSet,dimLocal,double,double > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;

      LagrangeFiniteElementMap ( const GridPart &gridPart, unsigned int order )
        : gridPart_( gridPart ), order_(order),
          localFeVector_( size(), nullptr ) {}

      static std::size_t size () { return LocalGeometryTypeIndex::size(dimLocal); }

      int order () const { return order_; }

      template< class Entity >
      int order ( const Entity &entity ) const { return order(); }

      template< class Entity >
      std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & >
      operator() ( const Entity &e ) const
      {
        unsigned int index = localFiniteElement(e.type());
        const LocalFiniteElementType &lfe = *(localFeVector_[index]);
        return std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & >
          { index, lfe.localBasis(), lfe.localInterpolation() };
      }

      bool hasCoefficients ( const GeometryType &type ) const { return true; }

      const LocalCoefficientsType& localCoefficients ( const GeometryType &type ) const
      {
        unsigned int index = localFiniteElement(type);
        return localFeVector_[index]->localCoefficients();
      }

      const GridPartType &gridPart () const { return gridPart_; }

    private:
      std::size_t localFiniteElement(const GeometryType &type) const
      {
        std::size_t index = LocalGeometryTypeIndex::index(type);
        if ( !localFeVector_[index] )
          localFeVector_[index] = new LocalFiniteElementType(type,order_);
        return index;
      }

      mutable std::vector<LocalFiniteElementType*> localFeVector_;
      const GridPartType &gridPart_;
      unsigned int order_;
    };

    template< class FunctionSpace, class GridPart,
              template< class, unsigned int > class PointSet = EquidistantPointSet,
              template< class > class Storage = CachingStorage >
    using LagrangeSpace
    = LocalFiniteElementSpace<
        LagrangeFiniteElementMap< FunctionSpace, GridPart, PointSet >,
        FunctionSpace, Storage >;
  }
}

#endif // HAVE_DUNE_LOCALFUNCTIONS
#endif // HAVE_DUNE_FEM_SPACE_LAGRANGE
