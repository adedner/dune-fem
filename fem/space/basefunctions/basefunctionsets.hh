#ifndef DUNE_BASEFUNCTIONSETS_HH
#define DUNE_BASEFUNCTIONSETS_HH

//- Dune includes 
#include <dune/common/fvector.hh>

//- local includes 
#include <dune/fem/space/common/basefunctioninterface.hh>
#include <dune/fem/space/common/basefunctionfactory.hh>
#include <dune/fem/space/common/dofstorage.hh>

namespace Dune {
/** @defgroup VectorialBaseFunction Vectorial Base Function Set
    @ingroup BaseFunction
@{
**/
  // Forward declarations
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet;
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VectorialBaseFunctionSet;

  //! Traits class for standard base function set
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct StandardBaseFunctionSetTraits
  {
    //! Export function space type
    typedef FunctionSpaceImp FunctionSpaceType;
    //! Type of the base function storage policy
    typedef StorageImp<FunctionSpaceType> StorageType;
    //! Exact type of the base function
    typedef StandardBaseFunctionSet<FunctionSpaceType, 
                                    StorageImp> BaseFunctionSetType;
    //! Factory type for the corresponding base functions (polymorphic)
    typedef BaseFunctionFactory<FunctionSpaceType> FactoryType;

  };

  //! \brief Standard base function set
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class StandardBaseFunctionSet : 
    public BaseFunctionSetDefault<StandardBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> >
  {
  public:
    typedef StandardBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> Traits;
    typedef typename FunctionSpaceImp::DomainType DomainType;
    typedef typename FunctionSpaceImp::RangeType RangeType;
    typedef typename FunctionSpaceImp::RangeFieldType DofType;
    typedef typename FunctionSpaceImp::JacobianRangeType JacobianRangeType;

  public:
    //! Constructor
    StandardBaseFunctionSet(const typename Traits::FactoryType& factory) :
      storage_(factory),
      diffVar0_(0),
      tmp_(0.),
      jTmp_(0.)
    {}

    //! Total number of base functions
    inline
    int numBaseFunctions() const;

    template <int diffOrd>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal,
                  RangeType& phi) const;


    //! This method should not be used directly. It is called by methods
    //! eval, jacobian, hessian, which appropriately set diffVariable.
    template <int diffOrd, class QuadratureType>
    inline
    void evaluate (int baseFunct, 
                   const FieldVector<int, diffOrd> &diffVariable, 
                   QuadratureType & quad, 
                   int quadPoint, RangeType & phi) const;

    template <class QuadratureType>
    inline
    DofType evaluateSingle(const int baseFunct,
                           const QuadratureType& quad, const int quadPoint,
                           const RangeType& factor) const;
      
    template <class Entity, class QuadratureType>
    inline
    DofType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const QuadratureType& quad, 
                                   const int quadPoint,
                                   const JacobianRangeType& factor) const;
  private:
    typename Traits::StorageType storage_;
    
    mutable FieldVector<int, 0> diffVar0_;
    mutable RangeType tmp_;
    mutable JacobianRangeType jTmp_;
  };

  //- VectorialBaseFunctionSet
  template <class FunctionSpaceImp, template <class> class StorageImp>
  struct VectorialBaseFunctionSetTraits 
  {
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef VectorialBaseFunctionSet<
      FunctionSpaceType, StorageImp> BaseFunctionSetType;
  };

  //! \brief Special base function implementation that takes advantage
  //! of the vectorial structure of the base functions.
  //! This base function can be used in conjunction with scalar basefunctions
  //! \f$ \phi_i \f$ which are extended to vectorial base functions like 
  //! \f$ \Phi_j = \phi_i e_k \f$, where \f$ e_k = [ \kronecker_ik ]_i \f$.
  template <class FunctionSpaceImp, template <class> class StorageImp>
  class VectorialBaseFunctionSet : 
    public BaseFunctionSetDefault<VectorialBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> >
  {
  private:
    typedef BaseFunctionSetDefault<
      VectorialBaseFunctionSetTraits<FunctionSpaceImp, StorageImp> > BaseType;
    typedef FunctionSpaceImp FunctionSpaceType;
    typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;
    typedef typename ScalarFunctionSpaceType::RangeType ScalarRangeType;
    typedef typename ScalarFunctionSpaceType::JacobianRangeType 
      ScalarJacobianRangeType;
    typedef BaseFunctionFactory<ScalarFunctionSpaceType> FactoryType;
    typedef typename FactoryType::BaseFunctionType BaseFunctionType;
    typedef StorageImp<ScalarFunctionSpaceType> StorageType;
    typedef VectorialBaseFunctionSetTraits<FunctionSpaceImp,StorageImp> Traits;
 
  public:
    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    typedef typename FunctionSpaceType::RangeFieldType DofType;
        
  public:
    //! Constructor
    VectorialBaseFunctionSet(const FactoryType& factory) :
      storage_(factory),
      util_(FunctionSpaceType::DimRange),
      tmp_(0),
      jTmp_(0) // changed to integer in case of integer func-space
//      tmp_(0.),
//      jTmp_(0.)
    {}

    ~VectorialBaseFunctionSet() {}

    inline
    int numBaseFunctions() const;

    template <int diffOrd>
    inline
    void evaluate(int baseFunct,
                  const FieldVector<int, diffOrd>& diffVar,
                  const DomainType& xLocal,
                  RangeType& phi) const;

    template <int diffOrd, class QuadratureType>
    inline
    void evaluate(int baseFunct, 
                  const FieldVector<deriType, diffOrd> &diffVariable, 
                  QuadratureType & quad, 
                  int quadPoint, RangeType & phi ) const;

    inline
    void jacobian(int baseFunct, const DomainType& xLocal, 
                  JacobianRangeType& gradPhi) const;

    template <class QuadratureImp>
    inline
    void jacobian(int baseFunct, QuadratureImp& quad, int quadPoint,
                  JacobianRangeType& gradPhi) const;

    inline
    DofType evaluateSingle(const int baseFunct, 
                           const DomainType& xLocal,
                           const RangeType& factor) const;
    
    template <class QuadratureType>
    inline
    DofType evaluateSingle(const int baseFunct,
                           const QuadratureType& quad, int quadPoint,
                           const RangeType& factor) const;
      
    template <class Entity>
    inline
    DofType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const DomainType& xLocal,
                                   const JacobianRangeType& factor) const;
    
    template <class Entity, class QuadratureType>
    inline
    DofType evaluateGradientSingle(const int baseFunct,
                                   const Entity& en,
                                   const QuadratureType& quad, int quadPoint,
                                   const JacobianRangeType& factor) const;

  private:
    StorageType storage_;
    DofConversionUtility<PointBased> util_;

    mutable FieldVector<int, 0> diffVar0_;
    mutable FieldVector<int, 1> diffVar1_;
    mutable ScalarRangeType tmp_;
    mutable ScalarJacobianRangeType jTmp_;
  };

/** @} **/
} // end namespace Dune

#include "basefunctionsets.cc"
#endif
