#ifndef DUNE_FEM_SPACE_FINITEVOLUME_DECLARATION_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_DECLARATION_HH

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeSpace
    // -----------------

    template< class FunctionSpace, class GridPart, int codim, class Storage >
    class FiniteVolumeSpace;

    // VertexCenteredFiniteVolumeSpace
    // -------------------------------

    template< class FunctionSpace, class GridPart, int codim, class Storage >
    class VertexCenteredFiniteVolumeSpace;



  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_DECLARATION_HH
