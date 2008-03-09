//************************************************************
//
//  (C) written and directed by Robert Kloefkorn 
//
//************************************************************
#include <iostream>
#include <vector>
#include <cassert>
#include <string>

#if HAVE_MPI == 1 
#error "Visualization only works without MPI" 
#endif 

#include <dune/common/misc.hh>
#include <dune/common/exceptions.hh>
using namespace Dune;

// include definition of grid type 
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

// include grape visualization 
#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/grid/io/visual/combinedgrapedisplay.hh>

// include data reading 
#include <dune/fem/io/visual/grape/datadisp/printhelp.cc>
#include <dune/fem/io/visual/grape/datadisp/readiotupledata.cc>
#include <dune/fem/io/visual/grape/datadisp/readioparams.cc> 
#include <dune/fem/io/parameter.hh>
#include <dune/fem/function/common/discretefunctionadapter.hh>

int main(int argc, char **argv)
{
  try {			         
    Parameter::append(argc,argv);
    if (argc < 2)
    {
      print_help(argv[0]);
      return(0);
    }   

    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "-help"))
    {
      print_help(argv[0]);
      return(0);
    }
    return readParameterList(argc,argv);
  }
  catch (Dune::Exception& e)
  {
    std::cerr << e << std::endl;
    return 1;
  }
  return 0;
}
