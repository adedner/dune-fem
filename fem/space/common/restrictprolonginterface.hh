#ifndef RESTRICTPROLONGINTERFACE
#define RESTRICTPROLONGINTERFACE

#include <dune/fem/quadrature/elementquadrature.hh>
#include <dune/fem/misc/combineinterface.hh>
namespace Dune{
/** @defgroup RestrictProlongInterface Restrict Prolong Interface
    @ingroup Adaptation
    Interface for restrict prolong operation on single elements.
  @{
 */

/*! \brief Interface class defining the local behaviour of the
    restrict/prolong operation (using BN)
 */
template <class RestProlImpTraits>
class RestrictProlongInterface {
public:  
  typedef typename RestProlImpTraits::RestProlImp RestProlImp;
  typedef typename RestProlImpTraits::RangeFieldType RangeFieldType;
  //! if weight is set, then ists assumend that we have always the same
  //! proportion between fahter and son volume 
  void setFatherChildWeight (const RangeFieldType& val) const {
    asImp().setFatherChildWeight(val);
  }
  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, 
		       bool initialize ) const {
    asImp().restrictLocal(father,son,initialize);
  }
  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, 
		      bool initialize ) const {
    asImp().prolongLocal(father,son,initialize);
  }

protected:  
  //! calculates the weight, i.e. (volume son)/(volume father)
  template <class EntityType>
  RangeFieldType calcWeight (EntityType &father, EntityType &son) const
  {
    assert( (son.geometry().volume() / father.geometry().volume()) > 0.0 );
    return (son.geometry().volume() / father.geometry().volume());
  }
  
private:
  RestProlImp& asImp() {
    return static_cast<RestProlImp&>(*this);
  }
  const RestProlImp& asImp() const {
    return static_cast<const RestProlImp&>(*this);
  }
};
template <class Base>
struct RestrictProlongTraits {
  typedef Base RestProlImp;
  typedef double RangeFieldType;
};
/*! \brief Allow the combination of two restrict/prolong instances
 */
template <class I1,class I2>
class RestrictProlongPair : 
    public RestrictProlongInterface<RestrictProlongTraits<RestrictProlongPair<I1,I2> > >,
    public PairOfInterfaces<I1,I2> {
public:  
  typedef PairOfInterfaces<I1,I2> BaseType;
  typedef typename BaseType::T1Type::RangeFieldType RangeFieldType;
  typedef typename BaseType::T2Type::RangeFieldType RangeFieldType2;

  RestrictProlongPair(I1 i1, I2 i2) : PairOfInterfaces<I1,I2>(i1,i2) {
    IsTrue<SameType<RangeFieldType,RangeFieldType2>::value>::yes();
  }
  
  //! if weight is set, then ists assumend that we have always the same
  //! proportion between fahter and son volume 
  void setFatherChildWeight (const RangeFieldType& val) const {
    this->first().setFatherChildWeight(val);
    this->second().setFatherChildWeight(val);    
  }
  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, 
		       bool initialize ) const {
    this->first().restrictLocal(father,son,initialize);
    this->second().restrictLocal(father,son,initialize);    
  }
  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, 
		      bool initialize ) const {
    this->first().prolongLocal(father,son,initialize);
    this->second().prolongLocal(father,son,initialize);    
  }
};
/** \brief This is a general restriction/prolongation operator
    which is speciallized for some discreteFunctionSpaces (e.g. DG)
 */
template <class DiscreteFunctionType> class RestrictProlongDefault :
    public RestrictProlongInterface<RestrictProlongTraits<RestrictProlongDefault<DiscreteFunctionType> > > {
};

/** \brief This is a simple restriction/prolongation operator for
 piecewise constant data stored on elements. 
*/
template <class DiscreteFunctionType>
 class RestProlOperatorFV : public RestrictProlongInterface<RestrictProlongTraits<RestProlOperatorFV<DiscreteFunctionType> > >
{
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType SpaceType; 
  typedef typename SpaceType :: GridPartType GridPartType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
public:  
  //! Constructor
  RestProlOperatorFV ( DiscreteFunctionType & df ) : df_ (df), weight_(-1.0)
  {
    // make sure index set can be used for adaptive computations 
  }

  //! if weight is set, then ists assumend that we have always the same
  //! proportion between fahter and son volume 
  void setFatherChildWeight (const RangeFieldType& val) const
  {
    // volume of son / volume of father  
    weight_ = val; 
  }
  
  //! restrict data to father 
  template <class EntityType>
  void restrictLocal ( EntityType &father, EntityType &son, bool initialize ) const
  {
    assert( df_.space().indexSet().adaptive() );

    assert( !father.isLeaf() );

    const RangeFieldType weight = (weight_ < 0.0) ? (this->calcWeight(father,son)) : weight_; 

    assert( weight > 0.0 );
    
    LocalFunctionType vati = df_.localFunction( father);
    LocalFunctionType sohn = df_.localFunction( son   );

    const int numDofs = vati.numDofs();
    if(initialize)
    {
      for(int i=0; i<numDofs; ++i)
      {
        vati[i] = weight * sohn[i];
      }
    }
    else 
    {
      for(int i=0; i<numDofs; ++i)
      {
        vati[i] += weight * sohn[i];
      }
    }
  }

  //! prolong data to children 
  template <class EntityType>
  void prolongLocal ( EntityType &father, EntityType &son, bool initialize ) const
  {
    assert( df_.space().indexSet().adaptive() );
    
    LocalFunctionType vati = df_.localFunction( father);
    LocalFunctionType sohn = df_.localFunction( son   );
    const int numDofs = vati.numDofs();
    for(int i=0; i<numDofs; ++i)
    {
      sohn[i] = vati[i];
    }
  }

private:
  mutable DiscreteFunctionType & df_;
  mutable RangeFieldType weight_;
};


/** @} end documentation group */

}
#endif
