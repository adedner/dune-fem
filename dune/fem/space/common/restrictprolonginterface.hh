#ifndef DUNE_RESTRICTPROLONGINTERFACE_HH
#define DUNE_RESTRICTPROLONGINTERFACE_HH

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/capabilities.hh>

//- local includes 
#include <dune/fem/misc/combineinterface.hh>
#include <dune/fem/gridpart/emptyindexset.hh>

namespace Dune{
/** @addtogroup RestrictProlongInterface 

    Interface for restriction and prolongation operation of data 
    on single elements.

    \remarks The Interface for a restriction and prolongation operation 
    is defined by the class RestrictProlongInterface.

    
  @{
 */

/*! @ingroup RestrictProlongInterface
    \brief Interface class defining the local behaviour of the
    restrict/prolong operation (using BN)

    \interfaceclass
 */
template <class RestProlImpTraits>
class RestrictProlongInterface {
public:  
  //! \brief type of restrict-prolong operator implementation 
  typedef typename RestProlImpTraits::RestProlImp RestProlImp;

  //! \brief field type of domain vector space
  typedef typename RestProlImpTraits::DomainFieldType DomainFieldType;

  /** \brief if weight is set, then its assumend 
      that we have always the same proportion between fahter and son volume 
      \param[in] weight proportion between fahter and son volume 
  */
  void setFatherChildWeight (const DomainFieldType& weight) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().setFatherChildWeight(weight));
  }
  
  /** \brief restrict data to father 
      \param[in] father Father Entity 
      \param[in] son Son Entity 
      \param[in] initialize <b>true</b> if restrictLocal is called for first time for this father
  */
  template <class EntityType>
  void restrictLocal ( EntityType &father, 
                       EntityType &son, 
            		       bool initialize ) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().restrictLocal(father,son,initialize));
  }
  
  /** \brief prolong data to children 
      \param[in] father Father Entity 
      \param[in] son Son Entity 
      \param[in] initialize <b>true</b> if restrictLocal is called for first time for this father
  */
  template <class EntityType>
  void prolongLocal ( EntityType &father, 
                      EntityType &son, 
            		      bool initialize ) const 
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().prolongLocal(father,son,initialize));
  }

  /** \brief add discrete function to communicator 
      \param[in] communicator Communicator to which internal discrete functions are added to 
  */
  template <class CommunicationManagerImp>
  void addToList(CommunicationManagerImp& communicator)
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().addToList(communicator));  
  }
  
protected:  
  /** \brief calculates the weight, i.e. (volume son)/(volume father)
      \param[in] father Father Entity 
      \param[in] son Son Entity 
      \return proportion between fahter and son volume
  */
  template< class EntityType >
  DomainFieldType calcWeight ( const EntityType &father, const EntityType &son ) const
  {
    const DomainFieldType weight = son.geometry().volume() / father.geometry().volume();
    assert( weight > DomainFieldType( 0 ) );
    return weight;
  }
  
private:
  RestProlImp& asImp() {
    return static_cast<RestProlImp&>(*this);
  }

  const RestProlImp& asImp() const {
    return static_cast<const RestProlImp&>(*this);
  }
};

/** \brief Traits class for derivation from RestrictProlongInterface. */
template <class Base, class DomainField>
struct RestrictProlongTraits {
  typedef Base RestProlImp;
  typedef DomainField DomainFieldType;
};

/*! \brief Allow the combination of two restrict/prolong instances
 */
template <class I1,class I2>
class RestrictProlongPair
: public RestrictProlongInterface< RestrictProlongTraits
    < RestrictProlongPair<I1,I2>,
      typename PairOfInterfaces< I1, I2 >::T1Type::DomainFieldType > >,
  public PairOfInterfaces<I1,I2>
{
  typedef PairOfInterfaces<I1,I2> BaseType;
public:  
  typedef typename BaseType::T1Type::DomainFieldType DomainFieldType;
  dune_static_assert( (Conversion< DomainFieldType, typename BaseType::T2Type::DomainFieldType >::sameType),
                      "DomainFieldType doesn't match." );

  RestrictProlongPair(I1 i1, I2 i2)
  : PairOfInterfaces<I1,I2>(i1,i2)
  {}
  
  //! if weight is set, then ists assumend that we have always the same
  //! proportion between fahter and son volume 
  void setFatherChildWeight (const DomainFieldType& val) const {
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
  
  //! prolong data to children 
  template <class CommunicatorImp>
  void addToList(CommunicatorImp& comm) 
  {
    this->first().addToList(comm); 
    this->second().addToList(comm);    
  }
};


/** \brief Interface default implementation for derived classes. 
    Note the difference to RestrictProlongDefaultImplementation, which
    represents the implementation for certain spaces. 
*/
template< class TraitsImp >
class RestrictProlongInterfaceDefault 
: public RestrictProlongInterface< TraitsImp >
{
protected:
  //! return true if father and son have the same index
  template< class IndexSetType, class EntityType >
  bool entitiesAreCopies ( const IndexSetType &indexSet,
                           const EntityType &father,
                           const EntityType &son ) const
  {
    assert( indexSet.persistent() );
    return (indexSet.index( father ) == indexSet.index( son ));
  }


public:
  typedef typename TraitsImp::DomainFieldType DomainFieldType;

  /** \brief explicit set volume ratio of son and father
   *
   *  \param[in]  weight  volume of son / volume of father
   *
   *  \note If this ratio is set, it is assume to be constant.
   */
  void setFatherChildWeight ( const DomainFieldType &weight ) const
  {
    // we do not use this information
  }
};


/** \brief This is a general restriction/prolongation operator
    which is speciallized for some discreteFunctionSpaces (e.g. DG)
    ,i.e., over the second template argument.
    The RestrictProlongDefault class is then one which should be
    applied by the user.
 */
template <class DiscreteFunctionType,class DiscreteFunctionSpace>
class RestrictProlongDefaultImplementation;

/** \brief This is a wrapper for the default implemented
    restriction/prolongation operator, which only takes a discrete
    function template 
 */
template <class DiscreteFunctionType> 
class RestrictProlongDefault
: public RestrictProlongDefaultImplementation<
           DiscreteFunctionType,
           typename DiscreteFunctionType::DiscreteFunctionSpaceType >
{
  typedef RestrictProlongDefaultImplementation<
           DiscreteFunctionType,
           typename DiscreteFunctionType::DiscreteFunctionSpaceType >
  BaseType;
public:
  RestrictProlongDefault(DiscreteFunctionType& discreteFunction) :
    BaseType(discreteFunction) {
    discreteFunction.enableDofCompression();
  }
private:
  RestrictProlongDefault();
};



/** \brief This is a simple restriction/prolongation operator for
 piecewise constant data stored on elements. 
*/
template< class DiscreteFunctionType >
class RestrictProlongPieceWiseConstantData
: public RestrictProlongInterfaceDefault
  < RestrictProlongTraits< RestrictProlongPieceWiseConstantData< DiscreteFunctionType >,
                           typename DiscreteFunctionType::DiscreteFunctionSpaceType::GridType::ctype > >
{
public:  
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType SpaceType; 
  typedef typename SpaceType :: GridPartType GridPartType;
  typedef typename SpaceType :: GridType GridType;

  typedef typename GridType::ctype DomainFieldType;

protected:
  typedef RestrictProlongInterfaceDefault
    < RestrictProlongTraits< RestrictProlongPieceWiseConstantData< DiscreteFunctionType >, DomainFieldType > >
    BaseType;

  using BaseType :: calcWeight;
  using BaseType :: entitiesAreCopies;

public:  
  //! Constructor
  explicit RestrictProlongPieceWiseConstantData( DiscreteFunctionType &df ) 
  : df_ (df), weight_( -1.0 )
  {}

  /** \brief explicit set volume ratio of son and father
   *
   *  \param[in]  weight  volume of son / volume of father
   *
   *  \note If this ratio is set, it is assume to be constant.
   */
  void setFatherChildWeight ( const DomainFieldType &weight ) const
  {
    weight_ = weight;
  }
  
  //! restrict data to father 
  template< class EntityType >
  void restrictLocal ( const EntityType &dad, const EntityType &filius, bool initialize ) const
  {
    // convert from grid entities to grid part entities 
    const EntityType& father = df_.space().gridPart().convert( dad );
    const EntityType& son    = df_.space().gridPart().convert( filius );

    // if father and son are copies, do nothing
    if( entitiesAreCopies( df_.space().indexSet(), father, son ) )
      return;

    assert( !father.isLeaf() );

    const DomainFieldType weight = (weight_ < 0.0) ? calcWeight( father, son ) : weight_; 

    assert( weight > 0.0 );
    
    LocalFunctionType lfFather = df_.localFunction( father );
    LocalFunctionType lfSon = df_.localFunction( son );

    const int numDofs = lfFather.numDofs();
    if( initialize )
    {
      for( int i = 0; i < numDofs; ++i )
        lfFather[ i ] = weight * lfSon[ i ];
    }
    else 
    {
      for( int i = 0; i < numDofs; ++i )
        lfFather[ i ] += weight * lfSon[ i ];
    }
  }

  //! prolong data to children 
  template< class EntityType >
  void prolongLocal ( const EntityType &dad, const EntityType &filius, bool initialize ) const
  {
    // convert from grid entities to grid part entities 
    const EntityType& father = df_.space().gridPart().convert( dad );
    const EntityType& son    = df_.space().gridPart().convert( filius );

    // if father and son are copies, do nothing
    if( entitiesAreCopies( df_.space().indexSet(), father, son ) )
      return;
    
    LocalFunctionType lfFather = df_.localFunction( father );
    LocalFunctionType lfSon = df_.localFunction( son );
    const int numDofs = lfFather.numDofs();
    for( int i = 0; i < numDofs; ++i )
      lfSon[ i ] = lfFather[ i ];
  } 

  //! add discrete function to communicator 
  template <class CommunicatorImp> 
  void addToList(CommunicatorImp& comm) 
  {
    comm.addToList(df_);
  }

private:
  mutable DiscreteFunctionType & df_;
  mutable DomainFieldType weight_;
};

/** \brief This is an empty restriction/prolongation operator 
*/
class RestrictProlongEmpty;
struct RestrictProlongEmptyTraits {
  typedef RestrictProlongEmpty RestProlImp;
  typedef double DomainFieldType;
};
class RestrictProlongEmpty
: public RestrictProlongInterfaceDefault
  < RestrictProlongEmptyTraits >
{
  typedef RestrictProlongInterfaceDefault
    < RestrictProlongEmptyTraits >
    BaseType;

public:  

protected:
  using BaseType :: calcWeight;
  using BaseType :: entitiesAreCopies;

public:  
  //! Constructor
  explicit RestrictProlongEmpty( ) 
  {}

  void setFatherChildWeight ( const double &weight ) const
  {}
  
  //! restrict data to father 
  template< class EntityType >
  void restrictLocal ( const EntityType &father, 
                       const EntityType &son, bool initialize ) const
  {}

  //! prolong data to children 
  template< class EntityType >
  void prolongLocal ( const EntityType &father, 
                      const EntityType &son, bool initialize ) const
  {}

  //! add discrete function to communicator 
  template <class CommunicatorImp> 
  void addToList(CommunicatorImp& comm) 
  {}

};
///@} 

} // end namespace Dune 

#endif
