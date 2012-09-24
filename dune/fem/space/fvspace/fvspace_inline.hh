#ifndef DUNE_FEM_FVSPACE_INLINE_HH
#define DUNE_FEM_FVSPACE_INLINE_HH

//- system includes 
#include <algorithm>

//- Dune includes 
#include <dune/fem/space/common/allgeomtypes.hh>

#include <dune/fem/space/fvspace.hh>

namespace Dune 
{

  namespace Fem 
  {

    template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
    inline FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
    FiniteVolumeSpace (GridPartType & gridPart,
            const InterfaceType commInterface,
            const CommunicationDirection commDirection)
    : DefaultType(gridPart, commInterface, commDirection),
      baseFuncSet_(),
      blockMapper_( BlockMapperProviderType::getObject( gridPart ) ),
      mapper_( blockMapper_ )
    {
      makeFunctionSpace(gridPart);
    }

    template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
    inline void FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
    makeFunctionSpace (GridPartType& gridPart) 
    {
      // search the macro grid for different element types 
      AllGeomTypes< typename GridPartType::IndexSetType ,
                    typename GridPartType::GridType > 
                allGeomTypes(gridPart.indexSet());
      
      const std::vector<GeometryType>& geomTypes = allGeomTypes.geomTypes(0);
      std::size_t maxDofs = 0;

      for( std::size_t i = 0; i < geomTypes.size(); ++i )
      {
        const GeometryType & geoType = geomTypes[i];
        if(baseFuncSet_.find( geoType ) == baseFuncSet_.end() ) 
        {
          // get pointer to base function set 
          const BaseFunctionSetImp* baseSet = 
            & SingletonProviderType::getObject(geoType);
          
          baseFuncSet_[ geoType ] = baseSet; 
          maxDofs = std::max( maxDofs, baseSet->size() );
        }
      }
      assert( std::size_t( blockMapper().maxNumDofs() * localBlockSize ) == maxDofs );
    }
      
    template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
    inline FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
    ~FiniteVolumeSpace () 
    {
      typedef typename BaseFunctionMapType :: iterator iterator;
      iterator end = baseFuncSet_.end();
      for(iterator it = baseFuncSet_.begin(); it != end; ++it)
      {
        const BaseFunctionSetImp * set = (*it).second; 
        SingletonProviderType::removeObject(*set);
      }

      baseFuncSet_.clear();
      BlockMapperProviderType::removeObject(blockMapper_);
    }  

    template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
    inline DFSpaceIdentifier FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
    type () const 
    {
      return FiniteVolumeSpace_id;
    }

    template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
    template <class EntityType> 
    inline const 
    typename FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::BaseFunctionSetType  
    FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
    baseFunctionSet (const EntityType &en) const 
    {
      return this->baseFunctionSet(en.type());
    }

    template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
    inline const 
    typename FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::BaseFunctionSetType  
    FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::
    baseFunctionSet (const GeometryType& geomType) const 
    {
      assert( baseFuncSet_.find(geomType) != baseFuncSet_.end());
      return BaseFunctionSetType(baseFuncSet_[geomType]);
    }

    template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
    inline typename FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::MapperType&
    FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::mapper() const {
      return mapper_;
    }
    template <class FunctionSpaceImp, class GridPartImp, int polOrd, template <class> class BaseFunctionStorageImp >
    inline typename FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::BlockMapperType&
    FiniteVolumeSpace<FunctionSpaceImp, GridPartImp, polOrd, BaseFunctionStorageImp>::blockMapper() const 
    {
      return blockMapper_;
    }
       
  } // namespace Fem  

} // namespace Dune 

#endif // #ifndef DUNE_FEM_FVSPACE_INLINE_HH
