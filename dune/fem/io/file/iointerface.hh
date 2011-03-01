#ifndef DUNE_IOINTERFACE_HH
#define DUNE_IOINTERFACE_HH

//- system includes 
#include <iostream>
#include <sys/stat.h>  
#include <sys/types.h>
#include <dirent.h>


//- Dune includes
#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

// defines function readParameter 
#include <dune/fem/io/file/asciiparser.hh>
// defines Parameter 
#include <dune/fem/io/parameter.hh>

// binary data io 
#include <dune/fem/io/io.hh>
#include <dune/fem/io/file/binarydataio.hh>

// input and output of tuples 
#include <dune/fem/io/file/iotuple.hh>

#include <dune/fem/misc/capabilities.hh>

// if grape was configured then include headers 
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

namespace Dune
{

/** @addtogroup DiscFuncIO
   The package dune-fem provides a number 
   of possibilities to write data to disk either
   for visualization and debugging purposes or
   for checkpointing. 
   Visualization output of discrete functions
   is provided at the moment for
   - GraPE   
   - VTK, e.g., for paraview
   - Matlab 
   .
   In addition data can be visualized in GraPE online during 
   a simulation by setting the parameter
   \b fem.io.grapedisplay to one.
   With the exception of Matlab output
   the Dune::DataWritter or the 
   Dune::Checkpointer 
   are used, the format is choosen through the parameter
   \b fem.io.outputformat;
   values are
   - 0: write data in GraPE format which can also
        be used for checkpointing - this format
        is basicly lossless.
   - 1: VTK cell data
   - 2: VTK vertex data
   .
   Utilities for Matlab output are availabe
   through the Dune::MatlabHelper.
  
   The \ref Dune::CheckPointer "checkpointing facility" 
   writes both the grid state, a number of
   discrete function, and other parameters
   provided by the user. This information can
   then be read from disc to continue the simulation
   from the saved state. Using the \ref Visualization
   "datadisp" programm the checkpoint files can also
   be used for visualization purposes.

   \remark The interface class for general output
           of DiscreteFunctions to disc is given
           by IOInterface. A general purpose data writter
           is provided by the class Dune::DataWriter.

   Since the GraPE output format is lossless
   it is also used by the Dune::CheckPointer
   class which also writes data to files
   but alternates between to filenames -
   while the Dune::DataWriter should be used
   to store data for postprocessing, the
   Dune::CheckPointer facility should be used 
   to be able to restart computations after
   a unexpected termination of a simulation.
   
   Data files are generated in the directory
   given by the \b fem.prefix parameter and
   with the file prefix chosen via the parameter
   \b fem.io.datafileprefix. The data to be
   saved to disk is given to the Dune::DataWriter instance
   through a reference to a Dune::Tuple of 
   discrete function pointer.
 
   For a time series, data can be either written for a fixed
   time step, or after a fixed number of iterations using the
   parameters
   \b fem.io.savestep or \b fem.io.savecount,
   respectivly.
   If a series of data is to be written without a real
   time variable available, e.g., a series of refined grids,
   startTime=0, endTime=100,
   \b fem.io.savestep=-1 and \b fem.io.savecount=1
   is a good choice to make; data is then written
   using \b datawriter.write(step,step).
 
   The following code snippet demonstrated the
   general usage of the Dune::DataWriter:
   \code
   typedef Dune::Tuple< DestinationType > IOTupleType;
   IOTupleType dataTup ( &U );
   typedef DataWriter< GridType, IOTupleType > DataWriterType;
   DataWriterType dataWriter( grid,gridfilename,dataTup,startTime,endTime );
   for (counter=0;time<endTime;counter++) {
     ...
     dataWriter.write(time,counter);
   }
   \endcode

    \femparam{fem.prefix, path used for all file output, ./}
    \femparam{fem.io.datafileprefix, prefix used for all data files}
    \femparam{fem.io.outputformat, output format, 0} 
              values are: 
                   - 0 = GRAPE (lossless format), 
                   - 1 = VTK, 
                   - 2 = VTK vertex data, 
                   - 3 = gnuplot
                   .
    \femparam{fem.io.grapedisplay, use grape for online visualization; default is 0 (no)}
    \femparam{fem.io.savestep, interval for writting data files}
     use value <0 to deativate
    \femparam{fem.io.savecount, number of time steps between writting
                                file}
    use value <0 to deactivate
**/


/** @ingroup DiscFuncIO
 \brief IOInterface to write data to hard disk 
 \interfaceclass
*/ 
class IOInterface {

protected: 
  //! default constructor 
  IOInterface() {}
  
public: 
  //! destructor 
  virtual ~IOInterface () {}

  //! return FEM key for macro grid reading 
  static std::string defaultGridKey( const int dimension , const bool  check = true )
  {
    const std::string oldGridKey( "fem.io.macroGridFile" );
    
    std::ostringstream gridKeyStream;
    gridKeyStream << oldGridKey << "_" << dimension << "d";
    const std::string newGridKey( gridKeyStream.str() );

    // check for old parameter 
    if( Parameter::exists( oldGridKey ) )
    {
      if( Parameter::exists( newGridKey ) )
      {
        std::cerr << "WARNING: ignoring `" << oldGridKey << "' because `"
                  << newGridKey << "' was also found in parameter file." << std::endl;
        return newGridKey;
      }
      else
      {
        std::cerr << "WARNING: change `" << oldGridKey << "' to `"  << newGridKey
                  << "' in parameter file." << std::endl;
        return oldGridKey;
      }
    }

    // check for parameter with dimension 
    if( check && !Parameter::exists( newGridKey ) )
    {
      std::cerr << "ERROR: Parameter `" << newGridKey << "' not found." << std::endl;
      DUNE_THROW( ParameterNotFound, "Parameter `" << newGridKey << "' not found." );
    }
    return newGridKey;
  }

  //! create given path in combination with rank 
  static void createPath ( const std::string &path ) DUNE_DEPRECATED
  {
    if( !createDirectory( path ) )
      std::cerr << "Failed to create path `" << path << "'." << std::endl;
  }
  
  //! create given path in combination with rank 
  static std::string createPathName(const std::string& pathPref, int rank )
  {
    std::string path(pathPref);
    
    // add proc number to path 
    {
      path += "_";
      std::stringstream rankDummy;
      rankDummy << rank; 
      path += rankDummy.str();
    }
    return path;
  }
  
  //! standard path reading and creation method 
  //! rank is added to output path 
  static std::string readPath()
  {
    return Parameter::commonOutputPath();
  }

  //! standard path reading and creation method 
  //! rank is added to output path 
  static std::string readPath(const std::string& paramfile) DUNE_DEPRECATED
  {
    std::string path;

    // read output path from parameter file 
    if( readParameter(paramfile,"OutputPath",path) ) 
    {
      return path;
    }
    else 
    {
      std::cerr << "Couldn't read output path, exiting... " << std::endl;
    }
    // default path is current directory 
    path = ".";
    return path;
  }
  
  /** \brief create global path for data output */
  template <class CommunicatorType>
  static void createGlobalPath(const CommunicatorType& comm,
          const std::string& path) 
  {
    // only rank 0 creates global dir 
    if( comm.rank() <= 0 )
    {
      // create directory 
      if( !createDirectory( path ) )
        std::cerr << "Failed to create path `" << path << "'." << std::endl;
    }

    // wait for all procs to arrive here 
    comm.barrier ();
  }

  // creates path and processor sub pathes 
  template <class CommunicatorType>
  static std::string createPath(const CommunicatorType& comm,
          const std::string& pathPrefix, 
          const std::string& dataPrefix,
          const int step)
  {
    // first proc creates directory 
    std::string filename(pathPrefix);
    filename += "/";
    filename += dataPrefix;
    std::string path = genFilename("",filename,step);

    // create global path 
    createGlobalPath( comm, path );
    
    // append path with p for proc 
    path += "/p";

    // create path if not exists 
    path = createPathName( path, comm.rank() );
    
    // create path if not exits 
    if( !createDirectory( path ) )
      std::cerr << "Failed to create path `" << path << "'." << std::endl;
    return path;
  }
  
  // creates path and processor sub pathes 
  static std::string createRecoverPath(
          const std::string& pathPrefix, 
          const int rank,
          const std::string& dataPrefix,
          const int step)
  {
    // first proc creates directory 
    std::string filename(pathPrefix);
    filename += "/";
    filename += dataPrefix;
    std::string path = genFilename("",filename,step);

    // append path with p for proc 
    path += "/p";

    // create proc dir 
    return createPathName( path , rank );
  }

  //! if grid is structured grid, write macro file 
  template <class GridImp>
  static void writeMacroGrid(const GridImp& grid, 
                             const std::string& macroname,
                             const std::string& path, 
                             const std::string& prefix) 
  {
    // do nothing for non-cartesian grids  
    if( ! Capabilities::isCartesian<GridImp>::v ) return;

    // create file descriptor 
    std::ifstream gridin(macroname.c_str());
    if( !gridin) 
    {
      std::cerr << "Couldn't open file `" << macroname << "' ! \n";
      return ;
    } 
        
    // read interval information of structured grid 
    dgf::IntervalBlock interval(gridin);
    if(!interval.isactive()) 
    {
      std::cerr<<"Did not find IntervalBlock in macro grid file `" << macroname << "' ! \n";
      return;
    }
    
    std::string filename(path);
    filename += "/";
    filename += prefix;
    filename += "_grid";

    saveCartesianGrid( grid, interval, filename); 
    return;
  }

  //! if grid is structured grid, write macro file 
  template <class GridImp>
  static void copyMacroGrid(const GridImp& g,
                            const std::string& orgPath,
                            const std::string& destPath, 
                            const std::string& prefix) 
  {
    // do nothing for unstructured grids 
    if( Capabilities::IsUnstructured<GridImp>::v ) return;

    {
      std::string filename(orgPath);
      filename += "/";
      filename += prefix;
      filename += "_grid";
      // add rank 
      filename += strRank(g.comm().rank());

      std::string destFilename(destPath);
      destFilename += "/";
      destFilename += prefix;
      destFilename += "_grid.macro";

      std::string cmd("cp ");
      cmd += filename; cmd += " ";
      cmd += destFilename;

      // copy file to actual path 
      if( system( cmd.c_str() ) < 0 )
      {
        std::cerr << "Unable to execute command: '" << cmd << "'." << std::endl;
        //DUNE_THROW( IOError, "Unable to execute command: '" << cmd << "'." );
      }
    }
  }

protected:
  //! create string containing rank 
  static std::string strRank(const int rank)
  {
    std::stringstream tmp;
    tmp << "." << rank;
    return tmp.str();
  }

  template <class Grid> 
  struct SaveParallelCartesianGrid
  {
    typedef FieldVector<int, Grid::dimension> iTupel;

    static void getCoordinates(const Grid&, iTupel& , iTupel&, iTupel& ,iTupel& )
    {
      DUNE_THROW(NotImplemented,"SaveParallelCartesianGrid not implemented for choosen GridType");
    }
  };

  template < int dim > 
  struct SaveParallelCartesianGrid< YaspGrid< dim > >
  {
    typedef YaspGrid< dim >  Grid;
    typedef FieldVector<int, Grid::dimension> iTupel;

    static void getCoordinates(const Grid& grid, const iTupel& anz, 
                               iTupel& origin, iTupel& originInterior, 
                               iTupel& lengthInterior )
    {
#if HAVE_MPI
      // Yasp only can do origin = 0
      origin = 0;

      enum { dimworld = Grid :: dimensionworld };
      enum { tag = MultiYGrid< dimworld, double> ::tag };
      YLoadBalance< dimworld > loadBalancer;
      Torus< dimworld > torus( MPI_COMM_WORLD, tag, anz, &loadBalancer );
      torus.partition( torus.rank(), origin, anz, originInterior, lengthInterior );
#endif
    }
  };

  template < class ct, int dim, SPRefinementStrategy strategy , class Comm > 
  struct SaveParallelCartesianGrid< SPGrid< ct, dim, strategy, Comm > >
  {
    typedef SPGrid< ct, dim, strategy, Comm > Grid;
    typedef FieldVector<int, Grid::dimension> iTupel;

    static void getCoordinates(const Grid& grid, const iTupel& anz, 
                               iTupel& origin, iTupel& originInterior, 
                               iTupel& lengthInterior )
    {
#if HAVE_MPI
      typedef SPMultiIndex< dim > SPGridMultiIndex ;
      SPGridMultiIndex begin = grid.gridLevel( 0 ).localMesh().begin();
      SPGridMultiIndex end   = grid.gridLevel( 0 ).localMesh().end();
      for( int i=0; i<dim; ++i) 
      {
        originInterior[ i ] = begin[ i ];
        lengthInterior[ i ] = end[ i ] - begin[ i ];
      }
#endif
    }
  };

  //! write my partition as macro grid 
  template <class GridImp>
  static void saveCartesianGrid (const GridImp& grid, 
                                 dgf::IntervalBlock& intervalBlock,
                                 std::string filename )
  {
    enum { dimworld = GridImp :: dimensionworld };
    const int rank = grid.comm().rank();

    FieldVector<double,dimworld> lang;
    FieldVector<int,dimworld>    anz;
    FieldVector<double,dimworld>  h;
    FieldVector<int,dimworld>    orig;
    
    if( intervalBlock.numIntervals() != 1 )
    {
      std::cerr << "Warning: Only 1 interval block is handled by "
                << "IOInterface::saveMacroGridImp" << std::endl;
    }

    typedef typename dgf::IntervalBlock::Interval Interval;
    const Interval &interval = intervalBlock.get( 0 );
    for( int i = 0; i < dimworld; ++i )
    {
      orig[ i ] = interval.p[ 0 ][ i ];
      lang[ i ] = interval.p[ 1 ][ i ] - interval.p[ 0 ][ i ];
      anz[ i ] = interval.n[ i ];
      h[ i ] = lang[ i ] / anz[ i ];
    }

    // write sub grid for this rank 
    {
      std::string subfilename (filename);
      // add rank 
      subfilename += strRank(rank);

#if HAVE_MPI 
      {
        typedef FieldVector<int,dimworld> iTupel;

        // origin is zero 
        iTupel o( orig );

        iTupel o_interior;
        iTupel s_interior;

        SaveParallelCartesianGrid< GridImp > :: 
          getCoordinates( grid, anz, o, o_interior, s_interior );

        FieldVector<double,dimworld> origin;
        for(int i=0; i<dimworld; ++i) 
          origin[ i ] = o[ i ];

        FieldVector<double,dimworld> sublang(0.0);
        for(int i=0; i<dimworld; ++i)
        {
          origin[i] = o_interior[i] * h[i];
          sublang[i] = origin[i] + (s_interior[i] * h[i]);
        }

        writeStructuredGrid(subfilename,origin,sublang,s_interior);
      }
#else
      {
        // in serial this should be zero 
        assert( rank == 0 );
        FieldVector<double,dimworld> zero(0.0);
        writeStructuredGrid(subfilename,zero,lang,anz);
      }
#endif
    }

    // write global grid on rank 0 
    if (rank == 0 )
    {
      // write global file for recovery 
      filename += ".global";
      FieldVector<double,dimworld> zero(0.0);
      writeStructuredGrid(filename,zero,lang,anz);
    }
  }

  //! write structured grid as DGF file 
  template <int dimworld>
  static void writeStructuredGrid(const std::string& filename,
                           const FieldVector<double,dimworld>& origin,
                           const FieldVector<double,dimworld>& lang,
                           const FieldVector<int,dimworld>& anz)
  {
    std::ofstream file (filename.c_str());
    if( file.is_open())
    {
      file << "DGF" << std::endl;
      file << "Interval" << std::endl;
      // write first point 
      for(int i=0;i<dimworld; ++i)
      {
        file << origin[i] << " ";
      }
      file << std::endl;
      // write second point 
      for(int i=0;i<dimworld; ++i)
      {
        file << lang[i] << " ";
      }
      file << std::endl;
      // write number of intervals in each direction 
      for(int i=0;i<dimworld; ++i)
      {
        file << anz[i] << " ";
      }
      file << std::endl;
      file << "#" << std::endl;

      file << "BoundaryDomain" << std::endl;
      file << "default 1" << std::endl;
      file << "#" << std::endl;
    }
    else
    {
      std::cerr << "Couldn't open file `" << filename << "' !\n";
    }
  }
}; // end class IOInterface 

} // end namespace Dune  
#endif
