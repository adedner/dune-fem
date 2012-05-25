#ifndef DUNE_FEM_ENTITYCOMMHELPER_HH
#define DUNE_FEM_ENTITYCOMMHELPER_HH

#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  template< InterfaceType iftype >
  struct EntityCommHelper;


  template<>
  struct EntityCommHelper< InteriorBorder_InteriorBorder_Interface >
  {
    static bool send ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity);
      return (p == BorderEntity);
    }

    static bool receive ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity);
      return (p == BorderEntity);
    }
  };


  template<>
  struct EntityCommHelper< InteriorBorder_All_Interface >
  {
    static bool send ( const PartitionType p )
    {
      return (p == InteriorEntity) || (p == BorderEntity);
    }

    static bool receive ( const PartitionType p )
    {
      //return true;
      return (p != InteriorEntity);
    }
  };


  template<>
  struct EntityCommHelper< Overlap_OverlapFront_Interface >
  {
    static bool send ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity) || (p == OverlapEntity);
      return (p != FrontEntity) && (p != GhostEntity);
    }

    static bool receive ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity) || (p == OverlapEntity) || (p == FrontEntity);
      return (p != GhostEntity);
    }
  };


  template<>
  struct EntityCommHelper< Overlap_All_Interface >
  {
    static bool send ( const PartitionType p )
    {
      //return (p == InteriorEntity) || (p == BorderEntity) || (p == OverlapEntity);
      return (p != FrontEntity) && (p != GhostEntity);
    }

    static bool receive ( const PartitionType p )
    {
      return true;
    }
  };


  template<>
  struct EntityCommHelper< All_All_Interface >
  {
    static bool send ( const PartitionType p )
    {
      return true;
    }

    static bool receive ( const PartitionType p )
    {
      return true;
    }
  };

}

#endif // #ifndef DUNE_FEM_ENTITYCOMMHELPER_HH
