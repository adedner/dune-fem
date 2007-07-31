#ifndef DUNE_DISCRETEFUNCTIONSPACE_HH
#define DUNE_DISCRETEFUNCTIONSPACE_HH

//- system includes
#include <assert.h>

//- Dune includes 
#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/fem/space/common/geometryconversion.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/basefunctioninterface.hh>


//- local includes 
#include "allgeomtypes.hh"
#include "singletonlist.hh"

namespace Dune{

  /** @defgroup DiscreteFunctionSpace DiscreteFunctionSpace
      @ingroup FunctionCommon
      This provides the interfaces for discrete function spaces. 
  
      @{
  */

  enum DFSpaceIdentifier {  LagrangeSpace_id , DGSpace_id , 
    CombinedSpace_id , FiniteVolumeSpace_id , DFAdapter_id };
 
  //**************************************************************************
  //
  //  --DiscreteFunctionSpaceInterface
  //
  /** \brief This is the interface for discrete function spaces. All methods
    declared here have to be implemented by the implementation class.
    The discrete function space always depends on a given grid. 
    For all diffrent element types of the grid the function space provides 
    a set of base functions for the different element types. 
    Because of the knowledge of on the one hand the grid an on the other
    hand the base functions sets, the discrete function space provides the size
    of the function space and a mapping from entity and local dof number
    to global dof number of the level of the entity.
    NOTE: A FunctionSpace is defined on a certain grid part.
  */
  template<class FunctionSpaceTraits>
  class DiscreteFunctionSpaceInterface : 
    public FunctionSpaceTraits::FunctionSpaceType  
  {
  public:
    //- Typedefs and enums
    //! type of traits class 
    typedef FunctionSpaceTraits Traits; 
    //! type of function space (define domain and range types)
    typedef typename FunctionSpaceTraits::FunctionSpaceType FunctionSpaceType;
    //! type of DiscretefunctionSapce implementation (Barton-Nackman)
    typedef typename FunctionSpaceTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    //! type of BaseFunctionSet of this space 
    typedef typename FunctionSpaceTraits::BaseFunctionSetType BaseFunctionSetType;
    //! type of underlying grid part 
    typedef typename FunctionSpaceTraits::GridPartType GridPartType;
    //! type of underlying grid  
    typedef typename GridPartType::GridType GridType;
    //! type of used index set 
    typedef typename GridPartType::IndexSetType IndexSetType;
    
    /** \brief iterator type traversing the set of 
        entities defining the discrete  function space 
        (only codim 0 at the moment, to be revised)
    */
    typedef typename GridPartType:: template Codim<0>::IteratorType IteratorType;
    
  public:
    //- Public methods
    //! Constructor 
    DiscreteFunctionSpaceInterface() : FunctionSpaceType() 
    {}

    /** \brief return type identifier of discrete function space 
        \return return type identifier of discrete function space
    */
    DFSpaceIdentifier type () const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().type());
      return asImp().type();
    }

    /** \brief get base function set for given entity 
        \param[in] entity Entity for which base function is requested 
        \return BaseFunctionSet for Entity 
    */
    template< class EntityType >
    const BaseFunctionSetType baseFunctionSet ( const EntityType &entity ) const 
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().baseFunctionSet( entity ));
      return asImp().baseFunctionSet( entity );
    }
  
    /** \brief return true if space contains global continuous functions 
       (i.e. for LagrangeSpace <b>true</b> is returned, for DiscontinuousGalerkinSpace <b>false</b> is returned. 
       \return <b>true</b> if space contians global continous functions, <b>false</b> otherwise 
    */
    bool continuous() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().continuous());
      return asImp().continuous(); 
    }

    /** \brief return reference to grid which belongs to discrete function space 
        \return reference to grid  
    */
    const GridType& grid() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().grid());
      return asImp().grid(); 
    }

    /** \brief return reference to grid which belongs to discrete function space 
        \return reference to grid  
    */
    GridType& grid() 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().grid());
      return asImp().grid(); 
    }

    /** \brief Return the corresponding grid part (const version) 
        \return reference to grid part 
    */ 
    const GridPartType& gridPart() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
      return asImp().gridPart(); 
    }

    /** \brief Return the corresponding grid part (const version) 
        \return reference to grid part 
    */ 
    GridPartType& gridPart() 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
      return asImp().gridPart(); 
    }

    /** \brief Return reference to the corresponding index set of the space 
        \return reference to index set  
    */ 
    const IndexSetType& indexSet() const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().indexSet());
      return asImp().indexSet(); 
    }

    /** \brief Return number of degrees of freedom for this space 
        \return number of degrees of freedom 
    */
    int size () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
      return asImp().size(); 
    }

    /** \brief For given entity map local dof number to global dof number 
        \param[in] entity Entity for which mapping is done 
        \param[in] localNum local dof number 
        \return global dof number, i.e. position in dof array 
    */    
    template <class EntityType>
    int mapToGlobal ( const EntityType &entity, 
                      const int localNum ) const
    {
      CHECK_INTERFACE_IMPLEMENTATION(asImp().mapToGlobal(entity, localNum ));
      return asImp().mapToGlobal ( entity , localNum );
    }

    /** \brief returns index of sequence in grid sequences 
        \return number of current sequence 
    */
    int sequence () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().sequence());
      return asImp().sequence(); 
    }

    /** \brief get global order of space  
        \return order of space, i.e. polynomial order of base functions 
    */
    int order () const 
    { 
      CHECK_INTERFACE_IMPLEMENTATION(asImp().order());
      return asImp().order(); 
    } 
  
    /** \brief Iterator pointing to first entity associated 
               with this discrete function space 
        \return Iterator pointing to first Entity
    */
    IteratorType begin() const {
      return asImp().begin();
    }

    /** \brief Iterator pointing behind last entity associated 
               with this discrete function space 
        \return Iterator pointing behind last Entity
    */
    IteratorType end() const {
      return asImp().end();
    }

  protected:
    //! Barton-Nackman trick 
    DiscreteFunctionSpaceType& asImp() 
    { 
      return static_cast<DiscreteFunctionSpaceType&>(*this); 
    }

    //! Barton-Nackman trick 
    const DiscreteFunctionSpaceType& asImp() const  
    { 
      return static_cast<const DiscreteFunctionSpaceType&>(*this); 
    }
  
  }; // end class DiscreteFunctionSpaceInterface

  //---------------------------------------------------------------------------
  //-
  //-  --DiscreteFunctionSpaceDefault
  //-
  //-
  //---------------------------------------------------------------------------
  /** \brief This is the class with default implementations for discrete function */
  template <class FunctionSpaceTraits>
  class DiscreteFunctionSpaceDefault :
    public DiscreteFunctionSpaceInterface<FunctionSpaceTraits>
  {
  public:
    typedef typename FunctionSpaceTraits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename FunctionSpaceTraits::GridPartType  GridPartType;

  public:
    //! Constructor
    DiscreteFunctionSpaceDefault(const GridPartType & gridPart) 
      : DiscreteFunctionSpaceInterface<FunctionSpaceTraits>() 
      , multipleGeometryTypes_( 
          AllGeomTypes< typename GridPartType::IndexSetType ,
                        typename GridPartType::GridType > :: multipleGeomTypes() )
    {
    }

    /** \brief returns true if grid has more than one geometry type (hybrid grid)
        \return <b>true</b> if  grid has more than one geometry type
        (hybrid grid), <b>false</b> otherwise 
    */
    bool multipleGeometryTypes() const { return multipleGeometryTypes_; }

    /** \brief returns true if base function sets depend on entity 
        \return <b>true</b> if base function set depend on entities, <b>false</b> otherwise 
    */
    bool multipleBaseFunctionSets() const { return false; }

  protected: 
    //! true if grid has more than one geometry type (hybrid grids)
    const bool multipleGeometryTypes_;
 
  private:
    //! Barton-Nackman trick 
    DiscreteFunctionSpaceType& asImp() 
    { 
      return static_cast<DiscreteFunctionSpaceType&>(*this); 
    }

    //! Barton-Nackman trick 
    const DiscreteFunctionSpaceType &asImp() const  
    { 
      return static_cast<const DiscreteFunctionSpaceType&>(*this); 
    }
  };

} // end namespace Dune 
#include <dune/fem/space/common/discretefunctionspacecommon.hh>
#endif
