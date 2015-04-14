#ifndef DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH
#define DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH

#include <dune/fem/misc/compatibility.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/lagrange/shapefunctionset.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>

#include "adaptmanager.hh"
#include "declaration.hh"
#include "generic.hh"
#include "mapper.hh"
#include "restrictprolong.hh"


namespace Dune
{

  namespace Fem
  {

    /** \addtogroup PAdaptiveLagrangeSpace
     *
     *  Provides access to base function sets for different element types in
     *  one grid and size of function space and maps from local to global dof
     *  number.
     *
     *  \note This space can only be used with special index sets. If you want
     *  to use the PAdaptiveLagrangeSpace with an index set only
     *  supporting the index set interface you will have to use the
     *  IndexSetWrapper class to provide the required functionality.
     *
     *  \note For adaptive calculations one has to use index sets that are
     *  capable of adaption (i.e. the method adaptive returns true). See also
     *  AdaptiveLeafIndexSet.
     */

    // PAdaptiveLagrangeSpaceTraits
    // ----------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    struct PAdaptiveLagrangeSpaceTraits
    {
      static_assert((polOrder > 0), "LagrangeSpace only defined for polOrder > 0" );

      typedef PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > DiscreteFunctionSpaceType;

      typedef FunctionSpace FunctionSpaceType;
      typedef GridPart GridPartType;

      static const int polynomialOrder = polOrder;

      static const bool continuousSpace = true ;
      static const int localBlockSize = FunctionSpaceType::dimRange;

      typedef PAdaptiveLagrangeMapper< GridPartType, polynomialOrder > BlockMapperType;

      typedef LagrangePointSet< GridPartType, polynomialOrder > CompiledLocalKeyType;

      static const int codimension = 0;

    private:
      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;

      static const int dimLocal = GridPartType::dimension;
      typedef typename FunctionSpace::ScalarFunctionSpaceType ScalarFunctionSpaceType;
      typedef typename ToNewDimDomainFunctionSpace< ScalarFunctionSpaceType, dimLocal >::Type ShapeFunctionSpaceType;

      typedef LagrangeShapeFunctionInterface< ShapeFunctionSpaceType > ShapeFunctionType;
      typedef SimpleShapeFunctionSet< ShapeFunctionType > SimpleShapeFunctionSetType;

    public:
      typedef SelectCachingShapeFunctionSet< SimpleShapeFunctionSetType, Storage > ScalarShapeFunctionSetType;

      template< int pOrd >
      struct ScalarShapeFunctionSetFactory
      {
        struct Type
        {
          static ScalarShapeFunctionSetType *createObject ( const GeometryType &type )
          {
            typedef LagrangeShapeFunctionFactory< ShapeFunctionSpaceType, polOrder > SimpleShapeFunctionSetFactoryType;
            return new ScalarShapeFunctionSetType( type, SimpleShapeFunctionSetType( SimpleShapeFunctionSetFactoryType( type ) ) );
          }

          static void deleteObject ( ScalarShapeFunctionSetType *object ) { delete object; }
        };
      };

      typedef ShapeFunctionSetProxy< ScalarShapeFunctionSetType > ScalarShapeFunctionSetProxyType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetProxyType, typename FunctionSpaceType::RangeType > ShapeFunctionSetType;

      typedef Dune::Fem::DefaultBasisFunctionSet< EntityType, ShapeFunctionSetType > BasisFunctionSetType;

      template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
        typedef Operation OperationType;
      };
    };



    // PAdaptiveLagrangeSpace
    // ----------------------

    /** \class   PAdaptiveLagrangeSpace
     *
     *  \ingroup PAdaptiveLagrangeSpace
     *
     *  \brief   Lagrange discrete function space
     */
    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    class PAdaptiveLagrangeSpace
    : public GenericDiscreteFunctionSpace< PAdaptiveLagrangeSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > >
    {
      typedef PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > ThisType;
      typedef GenericDiscreteFunctionSpace< PAdaptiveLagrangeSpaceTraits< FunctionSpace, GridPart, polOrder, Storage > > BaseType;

    public:
      typedef ThisType PAdaptiveLagrangeSpaceType;

      typedef typename BaseType::Traits Traits;

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename BaseType::IntersectionType IntersectionType;

      typedef typename BaseType::CompiledLocalKeyType CompiledLocalKeyType;
      typedef CompiledLocalKeyType LagrangePointSetType;

    protected:
      using BaseType::dfList_ ;
      using BaseType::searchFunction ;

    public:
      using BaseType::blockMapper;
      using BaseType::compiledLocalKey;
      using BaseType::continuous;
      using BaseType::gridPart;
      using BaseType::order;

      // default communication interface
      static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;
      // default communication direction
      static const CommunicationDirection defaultDirection = ForwardCommunication;

      /** \brief constructor
       *
       *  \param[in]  gridPart       grid part for the Lagrange space
       *  \param[in]  commInterface  communication interface to use (optional)
       *  \param[in]  commDirection  communication direction to use (optional)
       */
      explicit PAdaptiveLagrangeSpace ( GridPartType &gridPart,
                                        const InterfaceType commInterface = defaultInterface,
                                        const CommunicationDirection commDirection = defaultDirection )
      : BaseType( gridPart, commInterface, commDirection )
      {}

      // copy constructor needed for p-adaption
      PAdaptiveLagrangeSpace ( const PAdaptiveLagrangeSpace &other )
      : BaseType( other )
      {}

      /** @copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      bool continuous (const IntersectionType &intersection) const
      {
        if ( order() > 0 && intersection.conforming())
        {
          return true;
          if (intersection.neighbor())
            return (order(make_entity(intersection.inside())) == order(make_entity(intersection.outside())));
          else
            return true;
        }
        return false;
      }

      /** \brief provide access to the Lagrange point set for an entity
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the LagrangeDiscreteFunctionSpace.
       *
       *  \param[in]  entity  entity the Lagrange point set is requested for
       *
       *  \returns LagrangePointSet
       */
      template< class EntityType >
      const CompiledLocalKeyType &lagrangePointSet ( const EntityType &entity ) const
      {
        return compiledLocalKey( entity.type(),
                                 blockMapper().polynomOrder( entity ) );
      }

      /** \brief provide access to the Lagrange point set for a geometry type
       *
       *  \note This method is not part of the DiscreteFunctionSpaceInterface. It
       *        is unique to the LagrangeDiscreteFunctionSpace.
       *
       *  \param[in]  type  type of geometry the Lagrange point set is requested for
       *  \param[in]  order polynomial order the Lagrange point set is requested for
       *
       *  \returns LagrangePointSetType
       */
      const CompiledLocalKeyType &lagrangePointSet ( const GeometryType &type, const int order = BaseType::polynomialOrder ) const DUNE_DEPRECATED
      {
        return compiledLocalKey( type, order );
      }

      /** \brief add function to discrete function space for p-adaptation
       *         (currently only supported by AdaptiveDiscreteFunction )
       */
      template< class DiscreteFunction >
      void addFunction( DiscreteFunction &df ) const
      {
        assert( searchFunction( df ) == dfList_.end() );

        // select LagrangeInterpolation to be the LocalInterpolation
        typedef typename BaseType :: template PAdaptiveDiscreteFunctionEntry<
            DiscreteFunction, LagrangeInterpolation< DiscreteFunction, DiscreteFunction > > RealEntryType ;
        typedef typename BaseType :: PAdaptiveDiscreteFunctionEntryInterface
          EntryInterface;
        EntryInterface *entry = new RealEntryType( df );

        assert( entry );
        dfList_.push_front( entry );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVE_LAGRANGE_HH
