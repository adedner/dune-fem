#ifndef DUNE_FEM_PASSSTUB_HH
#define DUNE_FEM_PASSSTUB_HH

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/elementquadrature.hh>

#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/dgdiscretemodel.hh>
#include <dune/fem/quadrature/elementquadrature.hh>
namespace Dune
{

  namespace Fem
  {

    class ProblemStub;

    template< class GridImp >
    struct PassStubTraits
    {
      typedef GridImp GridType;
      typedef PassStubTraits<GridType> Traits;
      typedef ProblemStub DGDiscreteModelType;
      typedef FunctionSpace<double, double, GridType :: dimension, 1> FunctionSpaceType;
      typedef typename FunctionSpaceType :: RangeType  RangeType;
      typedef typename FunctionSpaceType :: DomainType DomainType;
      typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
      typedef LeafGridPart<GridType> GridPartType;
      typedef ElementQuadrature<GridPartType,0> VolumeQuadratureType;
      typedef ElementQuadrature<GridPartType,1> FaceQuadratureType;
      typedef LagrangeDiscreteFunctionSpace<
        FunctionSpaceType, GridPartType, 1> DiscreteFunctionSpaceType;
      typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    };

    template< class PreviousPass, int pId = -1 >
    class PassStub
    : public Fem::Pass< PassStubTraits< GridSelector::GridType >, PreviousPass, pId >
    {
      typedef PassStub< PreviousPass, pId > ThisType;
      typedef Fem::Pass< PassStubTraits< GridSelector::GridType >, PreviousPass, pId > BaseType;

    public:
      typedef PassStubTraits< GridSelector::GridType > PassStubTraitsType;
      typedef PreviousPass PreviousPassType;

      typedef typename BaseType::TotalArgumentType ArgumentType;
      typedef typename PassStubTraitsType::DestinationType DestinationType;

    public:
      explicit PassStub( PreviousPassType &previousPass )
      : BaseType( previousPass )
      {}
        
      virtual void compute(const ArgumentType& arg, DestinationType& dest) const
      {}

      virtual void allocateLocalMemory() {}
    };


    struct ProblemStub
    : public Fem::DGDiscreteModelDefaultWithInsideOutside< PassStubTraits< GridSelector::GridType >, 0 > 
    {
      typedef PassStubTraits<GridSelector::GridType> Traits;
    }; 

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASSSTUB_HH
