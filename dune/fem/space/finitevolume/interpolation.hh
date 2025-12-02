#ifndef DUNE_FEM_SPACE_FINITEVOLUME_INTERPOLATION_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_INTERPOLATION_HH

#include <functional>

#include <dune/fem/function/localfunction/average.hh>
#include <dune/fem/quadrature/cornerpointset.hh>

#include "basisfunctionset.hh"

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeLocalInterpolation
    // ------------------------------

    template< class GridPart, class Range >
    class FiniteVolumeLocalInterpolation
    {
      typedef FiniteVolumeLocalInterpolation< GridPart, Range > ThisType;

    public:
      /** \brief grid part type */
      typedef GridPart GridPartType;
      /** \brief entity type */
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      /** \brief basis function set type */
      typedef FiniteVolumeBasisFunctionSet< EntityType, Range > BasisFunctionSetType;

      /** \name Construction
       * \{
       */

      FiniteVolumeLocalInterpolation () {}

      void bind( const EntityType &entity ) {}
      void unbind() {}

      explicit FiniteVolumeLocalInterpolation ( const EntityType &entity )
      { bind(entity); }

      /** \} */

      /** \name Copying and assignment
       * \{
       */

      FiniteVolumeLocalInterpolation ( const ThisType & ) = default;

      FiniteVolumeLocalInterpolation &operator= ( const ThisType & ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \brief return basis function set */
      /*
      BasisFunctionSetType basisFunctionSet () const
      {
        return BasisFunctionSetType( entity() );
      }
      */

      /** \brief interpolate local function */
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        apply( localFunction, localDofVector );
      }

      /** \brief interpolate local function */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        // use local functions range type here in case those differ
        typename LocalFunction::RangeType value;
        LocalAverage< LocalFunction, GridPartType >::apply( localFunction, value );
        for( int i = 0; i < Range::dimension; ++i )
          localDofVector[ i ] = value[ i ];
      }

      /** \} */

    private:
      //const EntityType &entity () const { return entity_.get(); }

      //std::reference_wrapper< const EntityType > entity_;
    };



    // VertexCenteredFiniteVolumeLocalInterpolation
    // --------------------------------------------

    template< class DiscreteFunctionSpace >
    class VertexCenteredFiniteVolumeLocalInterpolation
    {
      typedef VertexCenteredFiniteVolumeLocalInterpolation< DiscreteFunctionSpace > ThisType;

      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

    public:
      /** \brief grid part type */
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      /** \brief entity type */
      typedef typename DiscreteFunctionSpaceType::EntityType EntityType;

      /** \brief type of range of space */
      typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

      /** \brief basis function set type */
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
      /** \brief basis function set type */
      typedef typename DiscreteFunctionSpaceType::BlockMapperType  BlockMapperType;

      /** \name Construction
       * \{
       */
      explicit VertexCenteredFiniteVolumeLocalInterpolation ( const DiscreteFunctionSpaceType &space )
        : space_( space )
      {}

      /** \} */

      /** \name Copying and assignment
       * \{
       */

      VertexCenteredFiniteVolumeLocalInterpolation ( const ThisType & ) = default;

      VertexCenteredFiniteVolumeLocalInterpolation &operator= ( const ThisType & ) = default;

      /** \} */

      /** \name Public member methods
       *  \{
       */

      void bind( const EntityType &entity ) {}
      void unbind() {}


      /** \brief return basis function set */
      /*
      BasisFunctionSetType basisFunctionSet () const
      {
        return BasisFunctionSetType( entity() );
      }
      */

      /** \brief interpolate local function */
      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        apply( localFunction, localDofVector );
      }

      /** \brief interpolate local function */
      template< class LocalFunction, class LocalDofVector >
      void apply ( const LocalFunction &localFunction, LocalDofVector &localDofVector ) const
      {
        typedef CornerPointSet< typename DiscreteFunctionSpaceType::GridPartType > CornerPointSetType;
        const auto& dualVolume = space_.dualVolumes();
        const auto& entity = localFunction.entity();

        const auto basisSet = space_.basisFunctionSet( entity );
        const int scalarSize = basisSet.scalarSize();

        const auto& blockMapper = space_.blockMapper();

        const auto vol = localFunction.geometry().volume() / double(scalarSize);

        std::vector< typename BlockMapperType::GlobalKeyType > entityDofs( scalarSize );
        blockMapper.map( entity, entityDofs );

        // use local functions range type here in case those differ
        typename LocalFunction::RangeType value;
        LocalAverage< LocalFunction, GridPartType >::apply( localFunction, value );

        CornerPointSetType quad( entity );
        assert( quad.nop() == scalarSize );

        for( int i=0, dof = 0; i<scalarSize; ++i )
        {
          localFunction.evaluate( quad.point( i ), value );//refElem.position( i, entity.dimension ), ;
          const double volFrac = vol / dualVolume[ entityDofs[ i ] ];
          for( int r = 0; r < RangeType::dimension; ++r, ++dof )
          {
            localDofVector[ dof ] += value[ r ] * volFrac ;
          }
        }
      }

      /** \} */

    protected:
      const DiscreteFunctionSpaceType& space_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_INTERPOLATION_HH
