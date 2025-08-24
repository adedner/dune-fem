#ifndef DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/function/localfunction/average.hh>
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/defaultcommhandler.hh>
#include <dune/fem/space/common/localrestrictprolong.hh>
#include <dune/fem/space/discontinuousgalerkin/generic.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>
#include <dune/fem/space/shapefunctionset/selectcaching.hh>

#include "basisfunctionsets.hh"
#include "dualbasisset.hh"
#include "declaration.hh"
#include "interpolation.hh"

namespace Dune
{

  namespace Fem
  {
    /////////////////////////////////////////////////////////
    //
    // Cell centered FV space
    //
    /////////////////////////////////////////////////////////


    // FiniteVolumeSpaceTraits
    // -----------------------

    template< class FunctionSpace, class GridPart, int codim, class Storage >
    struct FiniteVolumeSpaceTraits
    {
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > DiscreteFunctionSpaceType;

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static const int codimension = codim;

      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
      typedef FiniteVolumeBasisFunctionSets< EntityType, typename FunctionSpaceType::RangeType, FiniteVolumeBasisFunctionSet > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      typedef CodimensionMapper< GridPartType, codimension > BlockMapperType;
      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange > LocalBlockIndices;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Copy >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // FiniteVolumeSpace
    // -----------------

    template< class FunctionSpace, class GridPart, int codim = 0, class Storage = SimpleStorage >
    class FiniteVolumeSpace
    : public GenericDiscontinuousGalerkinSpace< FiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > >
    {
      typedef FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > ThisType;
      typedef GenericDiscontinuousGalerkinSpace< FiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > > BaseType;

    public:
      /** \brief maximum polynomial order of the space, here 0 since basis functions are constant */
      static const int polynomialOrder = 0;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::GridPartType */
      typedef typename BaseType::GridPartType GridPartType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::EntityType */
      typedef typename BaseType::EntityType EntityType;

      /** \brief basis function sets type */
      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      /** \brief local interpolation type */
      typedef FiniteVolumeLocalInterpolation< GridPart, typename BasisFunctionSetType::RangeType > InterpolationType;
      typedef InterpolationType LocalInterpolationType;
      typedef LocalInterpolationType InterpolationImplType;

      explicit FiniteVolumeSpace ( GridPartType &gridPart,
                                   const InterfaceType commInterface = InteriorBorder_All_Interface,
                                   const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, BasisFunctionSetsType(), commInterface, commDirection )
      {}

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      static DFSpaceIdentifier type () { return FiniteVolumeSpace_id; }

      /** \brief return local interpolation */
      InterpolationType interpolation () const
      {
        return InterpolationType();
      }

      /** \brief return local interpolation */
      static InterpolationType interpolation ( const EntityType &entity )
      {
        return InterpolationType();
      }

      LocalInterpolationType localInterpolation ( const EntityType &entity ) const
      {
        return LocalInterpolationType(entity);
      }

      /** \brief extend size of space beyond what the index set is delivering */
      void extendSize( const size_t extension ) { this->blockMapper().extendSize( extension ); }
    };


   /** \brief Local Mass Matrix for FV space */
    template <class FunctionSpaceImp, class GridPartImp, int polOrd,
              class BaseFunctionStorageImp,
              class VolumeQuadratureImp>
    class LocalMassMatrix<
      FiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >,
      VolumeQuadratureImp >
      : public LocalMassMatrixImplementationDgOrthoNormal<
          FiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >, VolumeQuadratureImp, false /* refElemScaling */>
    {
      typedef FiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceImp;
      typedef LocalMassMatrixImplementationDgOrthoNormal< DiscreteFunctionSpaceImp, VolumeQuadratureImp, false /* refElemScaling */ > BaseType;
    public:
      using BaseType :: BaseType;
    };





    // DefaultLocalRestrictProlong for FiniteVolumeSpace
    // -------------------------------------------------

    template< class FunctionSpace, class GridPart, int codim, class Storage >
    class DefaultLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    : public ConstantLocalRestrictProlong< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    {
    public:
      DefaultLocalRestrictProlong ( const FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > & )
      {}
    };



    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct hasFixedPolynomialOrder< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct hasStaticPolynomialOrder< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
        static const int order = 0;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isContinuous< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isHierarchic< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isLocalized< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false; // there is no method 'shapeFunctionSet( const EntityType & )'
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isAdaptive< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct threadSafe< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct viewThreadSafe< FiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities


    /////////////////////////////////////////////////////////
    //
    // Vertex centered FV space
    //
    /////////////////////////////////////////////////////////

    // VertexCenteredFiniteVolumeSpaceTraits
    // -------------------------------------

    template< class FunctionSpace, class GridPart, int codim, class Storage >
    struct VertexCenteredFiniteVolumeSpaceTraits
    {
      static_assert( codim == 0, "Only implemented for codim 0");

      typedef VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > DiscreteFunctionSpaceType;

      typedef GridPart GridPartType;
      typedef GridFunctionSpace< GridPartType, FunctionSpace > FunctionSpaceType;

      static const int codimension = 0;

      typedef typename GridPartType::template Codim< codimension >::EntityType EntityType;
      typedef FiniteVolumeBasisFunctionSets< EntityType, typename FunctionSpaceType::RangeType, VertexCenteredFiniteVolumeBasisFunctionSet > BasisFunctionSetsType;
      typedef typename BasisFunctionSetsType::BasisFunctionSetType BasisFunctionSetType;

      // select vertex unknowns (i.e. dimension)
      typedef CodimensionMapper< GridPartType, GridPartType::dimension > BlockMapperType;
      typedef Hybrid::IndexRange< int, FunctionSpaceType::dimRange > LocalBlockIndices;

      template <class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
      struct CommDataHandle
      {
        typedef Operation OperationType;
        typedef DefaultCommunicationHandler< DiscreteFunction, Operation > Type;
      };
    };



    // VertexCenteredFiniteVolumeSpace
    // -------------------------------

    template< class FunctionSpace, class GridPart, int codim = 0, class Storage = SimpleStorage >
    class VertexCenteredFiniteVolumeSpace
    : public GenericDiscontinuousGalerkinSpace< VertexCenteredFiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > >
    {
      typedef VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > ThisType;
      typedef GenericDiscontinuousGalerkinSpace< VertexCenteredFiniteVolumeSpaceTraits< FunctionSpace, GridPart, codim, Storage > > BaseType;

      typedef typename GridPart :: ctype ctype;

    public:
      using BaseType::begin;
      using BaseType::end;
      using BaseType::blockMapper;
      using BaseType::sequence;

      typedef std::vector< ctype > DualVolumesType;

      /** \brief maximum polynomial order of the space, here 0 since basis functions are constant */
      static const int polynomialOrder = 1;

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::GridPartType */
      typedef typename BaseType::GridPartType GridPartType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::EntityType */
      typedef typename BaseType::EntityType EntityType;

      /** \brief basis function sets type */
      typedef typename BaseType::BasisFunctionSetsType BasisFunctionSetsType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::BasisFunctionSetType */
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::BlockMapperType */
      typedef typename BaseType::BlockMapperType BlockMapperType;

      /** \brief local interpolation type */
      typedef VertexCenteredFiniteVolumeLocalInterpolation< ThisType > InterpolationType;
      typedef InterpolationType LocalInterpolationType;
      typedef LocalInterpolationType InterpolationImplType;

      explicit VertexCenteredFiniteVolumeSpace ( GridPartType &gridPart,
                                                 const InterfaceType commInterface = InteriorBorder_All_Interface,
                                                 const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, BasisFunctionSetsType(), commInterface, commDirection ),
          sequence_( -1 )
      {}

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::type */
      static DFSpaceIdentifier type () { return VertexCenteredFiniteVolumeSpace_id; }

      /** \copydoc Dune::Fem::DiscreteFunctionSpaceInterface::continuous */
      static constexpr bool continuous () { return true; }

      /** \brief return local interpolation */
      InterpolationType interpolation () const
      {
        return InterpolationType( *this );
      }

      /** \brief return local interpolation */
      [[deprecated]]
      static InterpolationType interpolation ( const EntityType &entity )
      {
        return interpolation ();
      }

      LocalInterpolationType localInterpolation ( const EntityType &entity ) const
      {
        return interpolation () ;
      }

      /** \brief extend size of space beyond what the index set is delivering */
      void extendSize( const size_t extension ) { this->blockMapper().extendSize( extension ); }

      const DualVolumesType& dualVolumes() const
      {
        const int spcSequence = sequence();
        if( sequence_ != spcSequence )
        {
          typedef typename BlockMapperType :: GlobalKeyType GlobalKeyType;
          dualVolumes_.clear();
          dualVolumes_.resize( blockMapper().size() );

          std::fill( dualVolumes_.begin(), dualVolumes_.end(), 0.0 );
          std::vector< GlobalKeyType > entityDofs( blockMapper().maxNumDofs() );

          const auto endit = end();
          for( auto it = begin(); it != endit; ++it )
          {
            const auto& entity = *it;
            const int numDofs = entity.subEntities( EntityType::dimension );

            entityDofs.resize( numDofs );
            blockMapper().map( entity, entityDofs );

            // this only works for triangles, fix for other cells
            const double dualVol = entity.geometry().volume() / double(numDofs);
            assert( entity.geometry().affine() );

            for( const auto& dualCell : entityDofs )
            {
              //std::cout << " dof " << dualCell << std::endl;
              dualVolumes_[ dualCell ] += dualVol;
            }
          }
          sequence_ = spcSequence;
        }
        return dualVolumes_;
      }

    protected:
      mutable DualVolumesType dualVolumes_;
      mutable int sequence_;
    };


#if 0
    // TODO, local mass inversion for dual cells
   /** \brief Local Mass Matrix for FV space */
    template <class FunctionSpaceImp, class GridPartImp, int polOrd,
              class BaseFunctionStorageImp,
              class VolumeQuadratureImp>
    class LocalMassMatrix<
      VertexCenteredFiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >,
      VolumeQuadratureImp >
      : public LocalMassMatrixImplementationDgOrthoNormal<
          VertexCenteredFiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp >, VolumeQuadratureImp, false /* refElemScaling */>
    {
      typedef VertexCenteredFiniteVolumeSpace< FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp > DiscreteFunctionSpaceImp;
      typedef LocalMassMatrixImplementationDgOrthoNormal< DiscreteFunctionSpaceImp, VolumeQuadratureImp, false /* refElemScaling */ > BaseType;
    public:
      using BaseType :: BaseType;
    };
#endif



    // DefaultLocalRestrictProlong for VertexCenteredFiniteVolumeSpace
    // -------------------------------------------------

#if 0
    // TODO, write dual cell restrict prolong operator
    template< class FunctionSpace, class GridPart, int codim, class Storage >
    class DefaultLocalRestrictProlong< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    : public ConstantLocalRestrictProlong< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
    {
    public:
      DefaultLocalRestrictProlong ( const VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > & )
      {}
    };
#endif




    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct hasFixedPolynomialOrder< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct hasStaticPolynomialOrder< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
        static const int order = 0;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isContinuous< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };

      /*
      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isHierarchic< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };
      */

      // TODO Default Quadrature

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isLocalized< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false; // there is no method 'shapeFunctionSet( const EntityType & )'
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct isAdaptive< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct threadSafe< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = false;
      };

      template< class FunctionSpace, class GridPart, int codim, class Storage >
      struct viewThreadSafe< VertexCenteredFiniteVolumeSpace< FunctionSpace, GridPart, codim, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities




  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_SPACE_HH
