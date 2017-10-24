#ifndef DUNE_FEM_ISTLLINEAROPERATOR_HH
#define DUNE_FEM_ISTLLINEAROPERATOR_HH

#if HAVE_DUNE_ISTL

// system includes
#include <string>

// local includes
#include <dune/fem/operator/matrix/istlmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    //! ISTLMatrixOperator
    template< class DomainFunction, class RangeFunction >
    struct ISTLLinearOperator
    : public ISTLMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
      typedef ISTLLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef ISTLMatrixObject< DomainSpaceType, RangeSpaceType > BaseType;

      typedef typename BaseType::LittleBlockType LittleBlockType;
      typedef typename LittleBlockType::field_type FieldType;

      static constexpr bool assembled = true;

      using BaseType::apply;
      using BaseType::communicate;
      using BaseType::matrix;

      //! constructor
      //! \param domainSpace space defining domain of operator
      //! \param rangeSpace  space defining range of operator
      //! \param param ISTL matrix parameters for preconditioning
      //!         - Preconditioning: {0,1,2,3,4,5,6} put -1 to get info
      //!         - Pre-iteration: number of iteration of preconditioner
      //!         - Pre-relaxation: relaxation factor
      ISTLLinearOperator( const std::string & ,
                          const DomainSpaceType &domainSpace,
                          const RangeSpaceType &rangeSpace,
                          const MatrixParameter& param = ISTLMatrixParameter() )
        : BaseType( domainSpace, rangeSpace, param )
      {}

      virtual void operator()( const DomainFunction &arg, RangeFunction &dest ) const
      {
        apply( arg, dest );
      }

      const BaseType &systemMatrix() const
      {
        return *this;
      }

      BaseType &systemMatrix()
      {
        return *this;
      }

      void maskRows ( const RangeFunction &maskFunction, FieldType diagonal = FieldType( 0 ) )
      {
        const auto &slaveDofs = maskFunction.space().slaveDofs();
        for( auto i = matrix().begin(), iend = matrix().end(); i != iend; ++i )
        {
          const auto &mask = maskFunction.dofVector()[ i.index() ];
          for( auto j = i->begin(), jend = i->end(); j != jend; ++j )
          {
            for( typename LittleBlockType::size_type k = 0; k < LittleBlockType::rows; ++k )
              (*j)[ k ] *= mask[ k ];

            if( (j.index() != i.index()) || slaveDofs.isSlave( i.index() ) )
              continue;

            for( typename LittleBlockType::size_type k = 0; k < LittleBlockType::rows; ++k )
              (*j)[ k ][ k ] += (FieldType( 1 ) - mask[ k ]) * diagonal;
          }
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLLINEAROPERATOR_HH
