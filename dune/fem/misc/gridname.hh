#ifndef DUNE_GRIDNAME_HH
#define DUNE_GRIDNAME_HH

#include <string>
#include <iostream>
#include <typeinfo>
#include <vector>
#include <climits>
#include <cstdlib>

namespace Dune 
{
  namespace Fem 
  {
    
    template <class GridImp> 
    struct GridName
    {
      static std::string computeName(const GridImp& grid) 
      {
        std::string name ( typeid( GridImp ).name() ); 

        size_t dunePos = name.find( "Dune" );
        name.erase( 0, dunePos+4 );

        char* endptr = 0;
        // get position of strings that are not numbers 
        long int result = strtol( name.c_str(), &endptr, 0 );
        if( result == LONG_MAX || result == LONG_MIN ) 
        {
          std::cerr << "GridName: faild to determine name of grid!" << std::endl;
          abort();
        }

        if( endptr ) 
          name = std::string( endptr );

        // Grid seems to be followed by IL 
        size_t pos = name.find( "GridI" );
        pos += 4; // add length of Grid to get pos of IL 

        if( pos < name.size() ) 
          name.erase( pos, name.size() - pos );
#ifndef NDEBUG 
        std::vector< std::string > knownGrids;
        knownGrids.push_back( "AlbertaGrid" );
        knownGrids.push_back( "ALUSimplexGrid" ); 
        knownGrids.push_back( "ALUConformGrid" );
        knownGrids.push_back( "ALUCubeGrid" ); 
        knownGrids.push_back( "ALUGrid" ); 
        knownGrids.push_back( "SGrid" );
        knownGrids.push_back( "YaspGrid" );
        knownGrids.push_back( "PrismGrid" );
        knownGrids.push_back( "UGGrid" );
        knownGrids.push_back( "OneDGrid" );
        knownGrids.push_back( "GeometryGrid" );
        knownGrids.push_back( "CartesianGrid" );
        knownGrids.push_back( "ParallelSimplexGrid" );
        knownGrids.push_back( "SPGrid" );
        knownGrids.push_back( "P4estGrid" );

        bool found = false ;
        for(size_t i=0; i<knownGrids.size(); ++i) 
        {
          if( name == knownGrids[ i ] )
          {
            found = true; 
            break;
          }
        }

        if( ! found ) 
        {
          std::cerr << "WARNING: Grid name `" << name << "' not found in list of known grids! Please add in file " << __FILE__ << std::endl;
        }
#endif

        return name;
      }

      static const std::string& name(const GridImp& grid) 
      {
        static std::string name ( computeName( grid ) ); 
        return name;
      }
    };

    template <class GridImp>
    inline std::string gridName(const GridImp& grid) 
    {
      return GridName< GridImp > :: name( grid ); 
    }
  }
}
#endif
