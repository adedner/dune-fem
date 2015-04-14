#ifndef DUNE_FEM_SPACE_PADAPTIVE_ADAPTMANAGER_HH
#define DUNE_FEM_SPACE_PADAPTIVE_ADAPTMANAGER_HH

#include <dune/grid/common/capabilities.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/localrestrictprolong.hh>

#include "declaration.hh"
#include "restrictprolong.hh"


namespace Dune
{

  namespace Fem
  {

    // DefaultLocalRestrictProlong
    // ---------------------------

    template< class FS, class GP, int ord, template< class > class S >
    class DefaultLocalRestrictProlong< Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >
    : public PLagrangeLocalRestrictProlong< typename GP::GridType, Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >
    {
    public:
      DefaultLocalRestrictProlong ( const Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > &space )
      : PLagrangeLocalRestrictProlong< typename GP::GridType, Fem::PAdaptiveLagrangeSpace< FS, GP, ord, S > >( space )
      {}
    };


    template< class FunctionSpaceImp, class GridPartImp, int polOrd, template< class > class StorageImp >
    class DefaultLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > >
    : public DiscontinuousGalerkinLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp >, false > // invert mass matrix or not
    {
    public:
      typedef DiscontinuousGalerkinLocalRestrictProlong< Fem::PAdaptiveDGSpace<
                                                         FunctionSpaceImp,
                                                         GridPartImp,
                                                         polOrd, StorageImp >, false > BaseType ;
      DefaultLocalRestrictProlong ( const Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, polOrd, StorageImp > & space )
        : BaseType( space )
      {}
    };

    template< class FunctionSpaceImp, class GridPartImp, template< class > class StorageImp >
    class DefaultLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    : public ConstantLocalRestrictProlong< Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > >
    {
    public:
      DefaultLocalRestrictProlong ( const Fem::PAdaptiveDGSpace< FunctionSpaceImp, GridPartImp, 0, StorageImp > & )
      {}
    };


    // pAdaptation
    // -----------

    template <class DF, class Vector, class DFS>
    void pAdaptation( DF& df, const Vector& polynomialOrders, const DFS &space, const int )
    {}

    /** \brief pAdaptation
        \param df  discrete function to adapt
        \param polynomialOrders  vector containing polynomial orders for each cell
        \param space  type of space tp be adapted
        \param polOrderShift possible shift of polynomial order (i.e. in case of
                             Taylor-Hood put -1 for the pressure) (default = 0)
    */
    template <class DF, class Vector,
              class FS, class GP, int p,
              template< class > class Storage >
    void pAdaptation( DF& df,
                      const Vector& polynomialOrders,
                      const Fem::PAdaptiveLagrangeSpace<FS,GP,p,Storage> &space,
                      const int polOrderShift = 0 )
    {
      /*
      typedef typename DF :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
      typedef typename GridPartType :: GridType  GridType;

      DiscreteFunctionSpaceType& newSpace = const_cast< DiscreteFunctionSpaceType& > (df.space());
      DiscreteFunctionSpaceType oldSpace( df.space().gridPart() );

      typedef DofManager< GridType > DofManagerType;
      DofManagerType& dm = DofManagerType :: instance( newSpace.grid() );

      for( const auto& entity : newSpace )
        oldSpace.blockMapper().setPolynomOrder( entity, newSpace.blockMapper().polynomOrder( entity ) );

      dm.resize();
      dm.compress();

      AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > tmp( "padaptation", oldSpace );

      newSpace.adapt( tmp )
      tmp.assign( df );

      for( const auto& entity : newSpace )
      {
        const int polOrder = polynomialOrders[ newSpace.indexSet().index( entity ) ] + polOrderShift ;
        newSpace.blockMapper().setPolynomOrder( entity, polOrder );
      }

      dm.resize();
      dm.compress();

      LagrangeInterpolation< DF > :: interpolateFunction( tmp, df );
      */
    }

    /** \brief pAdaptation
        \param df  discrete function to adapt
        \param polynomialOrders  vector containing polynomial orders for each cell
        \param polOrderShift possible shift of polynomial order (i.e. in case of
          Taylor-Hood put -1 for the pressure) (default = 0)
    */
    template <class DF, class Vector>
    void pAdaptation( DF& df,
                      const Vector& polynomialOrders,
                      const int polOrderShift = 0 )
    {
      pAdaptation( df, polynomialOrders, df.space(), polOrderShift );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_ADAPTMANAGER_HH
