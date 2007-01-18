#ifndef DUNE_DOFMAPPERINTERFACE_HH
#define DUNE_DOFMAPPERINTERFACE_HH

namespace Dune {

/** @defgroup DofMapper  DofMapperInterface 
 *  @ingroup DiscreteFunction

  @{
 */

//***********************************************************************
//
//  --MapperInterface 
//
//! Interface for calculating the size of a function space for a grid on a
//! specified level.
//! Furthermore the local to global mapping of dof number is done. 
//! Also during grid adaptation this mapper knows about old and new indices
//! of entities. 
//
//***********************************************************************
template <class DofMapperImp> 
class DofMapperInterface
{
public: 
  //! return number of dofs for special function space and grid on
  //! specified level
  int size () const 
  {
    return asImp().size();
  }

  //! map a local dof num of a given entity to a global dof num
  template <class EntityType>
  int mapToGlobal ( EntityType &en, int localNum ) const
  {
    return asImp().mapToGlobal( en , localNum );
  }

  //! return new size of space, i.e. after adaptation 
  int newSize() const 
  {
    return asImp().newSize();
  }
  
  //! return number of dofs on element
  int numDofs () const 
  {
    return asImp().numDofs();
  } 

  //! return number of holes in the data 
  int numberOfHoles() const 
  {
    return asImp().numberOfHoles(); 
  }
  
  //! return old index in dof array of given index ( for dof compress ) 
  int oldIndex (int num) const 
  { 
    return asImp().oldIndex(); 
  }
    
  //! return new index in dof array of given index ( for dof compress ) 
  int newIndex (int num) const 
  { 
    return asImp().newIndex(); 
  }

  //! return estimate for size additional need for restriction of data
  int additionalSizeEstimate() const 
  {  
    return asImp().additionalSizeEstimate(); 
  }

  //! return true if compress will affect data  
  bool needsCompress () const 
  {
    return asImp().needsCompress ();
  }
  
private:  
  //! Barton-Nackman trick 
  DofMapperImp &asImp()  { return static_cast<DofMapperImp &>(*this); };
  //! Barton-Nackman trick 
  const DofMapperImp &asImp() const { return static_cast<const DofMapperImp &>(*this); };
};

//! Default implementation for DofMappers, empty at this moment
template <class DofMapperImp> 
class DofMapperDefault : public DofMapperInterface<DofMapperImp>
{
  //! nothing here at the moment 
private:  
  //! Barton-Nackman trick 
  DofMapperImp &asImp()  { return static_cast<DofMapperImp &>(*this); };
  //! Barton-Nackman trick 
  const DofMapperImp &asImp() const { return static_cast<const DofMapperImp &>(*this); };
}; 

//! Key for Mapper singleton list 
template <class IndexSetImp>
class MapperSingletonKey 
{
  const IndexSetImp & indexSet_; 
  const int numDofs_; 
public:
  //! constructor taking index set and numDofs 
  MapperSingletonKey(const IndexSetImp & indexSet, int numDofs ) : indexSet_(indexSet) ,  numDofs_(numDofs) {}
  //! copy constructor 
  MapperSingletonKey(const MapperSingletonKey &org) : indexSet_(org.indexSet_) , numDofs_(org.numDofs_) {}
  //! returns true if indexSet pointer and numDofs are equal 
  bool operator == (const MapperSingletonKey & otherKey) const 
  {
    return ((&indexSet_ == &otherKey.indexSet_) && (numDofs_ == otherKey.numDofs_));
  }

  //! return reference to index set 
  const IndexSetImp & indexSet() const { return indexSet_; }
  //! return number of dofs 
  const int numDofs () const { return numDofs_; }
};

//! Factory class for SingletonList to tell how objects are created and
//! how compared.
template <class KeyImp, class ObjectImp>
class MapperSingletonFactory
{
  public:
  //! create new mapper  
  static ObjectImp * createObject( const KeyImp & key )
  {
    // create Object of MapperType = ObjectImp 
    return new ObjectImp(key.indexSet(),key.numDofs());
  }
  //! delete mapper object 
  static void deleteObject( ObjectImp * obj )
  {
    delete obj;
  }
};

/** @} end documentation group */

} // end namespace Dune
#endif
