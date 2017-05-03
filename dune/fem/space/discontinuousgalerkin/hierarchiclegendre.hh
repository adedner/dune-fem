#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH

#include <cassert>

#include <dune/common/power.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/legendre.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

#include "basisfunctionsets.hh"
#include "generic.hh"
#include "interpolation.hh"
#include "shapefunctionsets.hh"

namespace Dune
{

  namespace Fem
  {

    // HierarchicLegendreDiscontinuousGalerkinSpaceTraits
    // --------------------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct HierarchicLegendreDiscontinuousGalerkinSpaceTraits
    {
      typedef HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int codimension = 0;

      typedef Dune::Fem::FunctionSpace<
          typename FunctionSpace::DomainFieldType, typename FunctionSpace::RangeFieldType,
           GridPartType::dimension, 1
        > ScalarShapeFunctionSpaceType;

      struct ScalarShapeFunctionSet
        : public Dune::Fem::HierarchicLegendreShapeFunctionSet< ScalarShapeFunctionSpaceType >
      {
        typedef Dune::Fem::HierarchicLegendreShapeFunctionSet< ScalarShapeFunctionSpaceType > BaseType;

        static const int numberShapeFunctions =
            StaticPower<polOrder+1,ScalarShapeFunctionSpaceType::dimDomain>::power;
      public:
        explicit ScalarShapeFunctionSet ( Dune::GeometryType type )
          : BaseType( polOrder )
        {
          assert( type.isCube() );
          assert( size() == BaseType::size() );
        }

        // overload size method because it's a static value
        unsigned int size() const { return numberShapeFunctions; }
      };

      typedef SelectCachingShapeFunctionSets< GridPartType, ScalarShapeFunctionSet, Storage > ScalarShapeFunctionSetsType;
      typedef VectorialShapeFunctionSets< ScalarShapeFunctionSetsType, typename FunctionSpaceType::RangeType > ShapeFunctionSetsType;

      typedef DefaultBasisFunctionSets< GridPartType, ShapeFunctionSetsType > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;
      static const int localBlockSize
        = FunctionSpaceType::dimRange*StaticPower< polOrder+1, GridPartType::dimension >::power;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // HierarchicLegendreDiscontinuousGalerkinSpace
    // --------------------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class HierarchicLegendreDiscontinuousGalerkinSpace
    : public GenericDiscontinuousGalerkinSpace< HierarchicLegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef GenericDiscontinuousGalerkinSpace< HierarchicLegendreDiscontinuousGalerkinSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      using BaseType::basisFunctionSet;

      static const int polynomialOrder = polOrder;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      typedef DiscontinuousGalerkinLocalL2Projection< GridPartType, BasisFunctionSetType > InterpolationType;

      explicit HierarchicLegendreDiscontinuousGalerkinSpace ( GridPartType &gridPart,
                                                              const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                              const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, makeBasisFunctionSets( gridPart ), commInterface, commDirection )
      {}

      static DFSpaceIdentifier type () { return HierarchicLegendreDGSpace_id; }

      InterpolationType interpolation ( const EntityType &entity ) const
      {
        return InterpolationType( basisFunctionSet( entity ) );
      }

    private:
      static BasisFunctionSetsType makeBasisFunctionSets ( const GridPartType &gridPart )
      {
        typedef typename BasisFunctionSetsType::ShapeFunctionSetsType ShapeFunctionSetsType;
        ShapeFunctionSetsType shapeFunctionSets( gridPart );
        return BasisFunctionSetsType( std::move( shapeFunctionSets ) );
      }
    };



    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasFixedPolynomialOrder< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasStaticPolynomialOrder< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isContinuous< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isLocalized< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isAdaptive< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct threadSafe< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct viewThreadSafe< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isHierarchic< HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_HLEGENDRE_HH
