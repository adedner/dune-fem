#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_SELECTCACHING_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_SELECTCACHING_HH

// dune-fem includes
#include <dune/fem/space/shapefunctionset/caching.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class >
    class CachingStorage;

    template< class >
    class SimpleStorage;



    // SelectCachingShapeFunctionSet
    // -----------------------------

    template< class ShapeFunctionSet, template< class > class Storage >
    class SelectCachingShapeFunctionSet;

    template< class ShapeFunctionSet >
    class SelectCachingShapeFunctionSet< ShapeFunctionSet, CachingStorage >
    : public CachingShapeFunctionSet< ShapeFunctionSet >
    {
      typedef CachingShapeFunctionSet< ShapeFunctionSet > BaseType;

    public:
      explicit SelectCachingShapeFunctionSet ( const GeometryType &type,
                                               const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : BaseType( type, shapeFunctionSet )
      {}
    };

    template< class ShapeFunctionSet >
    class SelectCachingShapeFunctionSet< ShapeFunctionSet, SimpleStorage >
    : public ShapeFunctionSet
    {
      typedef ShapeFunctionSet BaseType;

    public:
      explicit SelectCachingShapeFunctionSet ( const GeometryType &type,
                                               const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : BaseType( shapeFunctionSet )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_SELECTCACHING_HH