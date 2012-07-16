#ifndef DUNE_GENERICBASEFUNCTIONSETS_HH
#define DUNE_GENERICBASEFUNCTIONSETS_HH

#include <vector>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>

namespace Dune
{
  
  namespace Fem 
  {

  // Internal Forward Declarations
  // -----------------------------

  template< class LocalBasis >
  class GenericBaseFunctionSet;



  // GenericBaseFunctionSetTraits
  // ----------------------------

  template< class LocalBasis >
  class GenericBaseFunctionSetTraits
  {
    typedef typename LocalBasis::Traits::DomainFieldType DomainFieldType;
    typedef typename LocalBasis::Traits::RangeFieldType RangeFieldType;

    static const int dimDomain = LocalBasis::Traits::dimDomain;
    static const int dimRange = LocalBasis::Traits::dimRange;

  public:
    typedef FunctionSpace< DomainFieldType, RangeFieldType, dimDomain, dimRange >
      FunctionSpaceType;
    typedef GenericBaseFunctionSet< LocalBasis > BaseFunctionSetType;
  };



  // GenericBaseFunctionSet
  // ----------------------

  template< class LocalBasis >
  class GenericBaseFunctionSet
  : public BaseFunctionSetDefault< GenericBaseFunctionSetTraits< LocalBasis > >
  {
    typedef GenericBaseFunctionSet< LocalBasis > ThisType;
    typedef BaseFunctionSetDefault< GenericBaseFunctionSetTraits< LocalBasis > > BaseType;

  public:
    typedef LocalBasis LocalBasisType;

    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;

  public:
    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType RangeType;

    static const int dimDomain = BaseType::dimDomain;
    static const int dimRange = BaseType::dimRange;

    typedef typename LocalBasisType::Traits::JacobianType JacobianRangeType;

  public:
    // use evaluate of default implementation 
    using BaseType::evaluate;
    using BaseType::jacobian;

  public:
    //! Constructor
    explicit GenericBaseFunctionSet ( const LocalBasis &localBasis,
                                      const GeometryType &geometryType )
    : localBasis_( localBasis ),
      geometryType_( geometryType )
    {}

    /** \copydoc Dune::Fem::BaseFunctionSetInterface::size */
    size_t size () const
    {
      return localBasis_.size();
    }
    
    /** \copydoc Dune::Fem::BaseFunctionSetInterface::geometryType */
    GeometryType geometryType () const
    {
      return geometryType_;
    }
 
    /** \copydoc Dune::Fem::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<int,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */ 
    template< int diffOrd, class PointType >
    void evaluate ( const int baseFunction,
                    const FieldVector< int, diffOrd > &diffVariable,
                    const PointType &x,
                    RangeType &phi ) const
    {
      std::cerr << "derivates of order >= 2 not supported." << std::endl;
      abort();
    }

    template< class PointType >
    void evaluate ( const int baseFunction,
                    const FieldVector< int, 0 > &diffVariable,
                    const PointType &x,
                    RangeType &phi ) const
    {
      static std::vector< RangeType > tmpphi;
      localBasis_.evaluateFunction( coordinate( x ), tmpphi );
      assert( (unsigned int)baseFunction < tmpphi.size() );
      phi = tmpphi[ baseFunction ];
    }

    template< class PointType >
    void evaluate ( const int baseFunction,
                    const FieldVector< int, 1 > &diffVariable,
                    const PointType &x,
                    RangeType &phi ) const
    {
      static std::vector< JacobianRangeType > tmpphi;
      localBasis_.evaluateJacobian( coordinate( x ), tmpphi );
      assert( (unsigned int)baseFunction < tmpphi.size() );
      for( int j = 0 ; j < dimRange ; ++j )
        phi[ j ] = tmpphi[ baseFunction ][ j ][ diffVariable[ 0 ] ];
    }

    /** \copydoc Dune::Fem::BaseFunctionSetInterface::jacobian(const int baseFunction,const PointType &x,JacobianRangeType &phi) const */ 
    template< class PointType >
    void jacobian ( const int baseFunction,
                    const PointType &x,
                    JacobianRangeType &phi ) const
    {
      static std::vector< JacobianRangeType > tmpphi;
      localBasis_.evaluateJacobian( coordinate( x ), tmpphi );
      assert( (unsigned int)baseFunction < tmpphi.size() );
      phi = tmpphi[ baseFunction ];
    }
    
  private:
    GenericBaseFunctionSet ( const ThisType & );

    const LocalBasisType &localBasis_;
    GeometryType geometryType_;
  };

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_GENERICBASEFUNCTIONSETS_HH
