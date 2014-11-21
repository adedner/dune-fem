#ifndef DUNE_FEM_SPOPERATOR_HH
#define DUNE_FEM_SPOPERATOR_HH

#include <dune/fem/operator/matrix/spmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    // SparseRowLinearOperator
    // -----------------------

    template< class DomainFunction, class RangeFunction >
    class SparseRowLinearOperator
    : public SparseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public Fem::AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef SparseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType > Base;

    public:
      typedef typename Base::DomainSpaceType DomainSpaceType;
      typedef typename Base::RangeSpaceType RangeSpaceType;

      /** \copydoc Fem::Operator::assembled */
      static const bool assembled = true ;

      using Base::apply;

      SparseRowLinearOperator ( const std::string &name,
                                const DomainSpaceType &domainSpace,
                                const RangeSpaceType &rangeSpace,
                                const std::string &paramfile = "" )
      : Base( domainSpace, rangeSpace, paramfile )
      {}

      virtual void operator() ( const DomainFunction &arg, RangeFunction &dest ) const
      {
        Base::apply( arg, dest );
      }

      const Base &systemMatrix () const
      {
        return *this;
      }

      void communicate () const
      {
      }
    };
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPLINEAR_HH
