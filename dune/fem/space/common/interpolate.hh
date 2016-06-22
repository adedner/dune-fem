#ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
#define DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <vector>

namespace Dune
{

  namespace Fem
  {

    // interpolate
    // -----------

    template< class GridFunction, class DiscreteFunction, unsigned int partitions >
    static inline void interpolate ( const GridFunction &u, DiscreteFunction &v, PartitionSet< partitions > ps )
    {
      // reserve memory for local dof vector
      std::vector< typename DiscreteFunction::RangeFieldType > ldv;
      ldv.reserve( v.space().blockMapper().maxNumDofs() * DiscreteFunction::DiscreteFunctionSpaceType::localBlockSize );

      typename GridFunction::LocalFunctionType uLocal( u );

      // iterate over selected partition
      for( const auto entity : elements( v.gridPart(), ps ) )
      {
        // obtain local interpolation
        const auto interpolation = v.space().interpolation( entity );

        // resize local dof vector
        ldv.resize( v.space().basisFunctionSet( entity ).size() );

        // interpolate u locally
        uLocal.init( entity );
        interpolation( uLocal, ldv );

        // write local dofs into v
        v.setLocalDofs( entity, ldv );
      }
    }

    /**
     * \function interpolate
     * \ingroup  DiscreteFunctionSpace
     * \brief    perform native interpolation of a discrete function space
     *
     * By definition of its degrees of freedom, each discrete function space
     * has a native interpolation, which can be computed very quickly.
     *
     * For example, the native interpolation of a Lagrange discrete function
     * space is the evaluation in its Lagrange points.
     * An orthonormal DG space would instead perform an \f$L^2\f$-Projection.
     *
     * The actual implementation must locally be provided by the discrete
     * function space through the method
     * \code
     * template< class LocalFunction, class LocalDofVector >
     * void interpolate ( const LocalFunction &f, LocalDofVector &dofs ) const;
     * \endcode
     *
     * \param[in]   u  grid function to interpolate
     * \param[out]  v  discrete function to represent the interpolation
     */
    template< class GridFunction, class DiscreteFunction >
    static inline void interpolate ( const GridFunction &u, DiscreteFunction &v )
    {
      // just call interpolate for the all partition
      interpolate( u, v, Partitions::all );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_COMMON_INTERPOLATE_HH
