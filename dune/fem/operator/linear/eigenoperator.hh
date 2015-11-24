#ifndef DUNE_FEM_EIGENOPERATOR_HH
#define DUNE_FEM_EIGENOPERATOR_HH

#ifdef HAVE_EIGEN

// system includes
#include <string>

// local includes
#include <dune/fem/operator/matrix/eigenmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    //! EigenLinearOperator
    template< class DomainFunction, class RangeFunction >
    struct EigenLinearOperator
    : public EigenMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public Fem::AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef typename DomainFunction::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType RangeSpaceType;
      typedef EigeninearOperator< DomainFunction, RangeFunction > ThisType;
      typedef EigenMatrixObject< DomainSpaceType, RangeSpaceType > BaseType;

      static constexpr bool assembled = true ;

      using BaseType::apply;
      using BaseType::communicate

      EigenLinearOperator( const std::string & ,
                           const DomainSpaceType &domainSpace,
                           const RangeSpaceType &rangeSpace,
                           const std::string &paramfile = "" ) :
        BaseType( domainSpace, rangeSpace, paramfile )
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
    };
  } // namespace Fem

} // namespace Dune

#endif

#endif // #ifndef DUNE_FEM_SPLINEAR_HH