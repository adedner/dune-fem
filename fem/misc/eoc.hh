#ifndef DUNE_EOCOUTPUT_HH
#define DUNE_EOCOUTPUT_HH

namespace Dune{

/**  
    @ingroup HelperClasses
    \brief Write a self contained tex table and a gnuplot
    file for eoc runs with timing information.
    
    Constructor takes base name (filename) of file and
    generates two files:
    filename.tex, filename_body.tex and filename.gnu
    The file filename_body.tex hold the actual body
    of the tex file which is included in filename.tex
    but can also be used to combine runs with different
    parameters.
 */
class EocOutput {
 std::string outputFile;
  int level;
  std::vector<double> prevError;
  double prevSize;
  bool initial;
  bool finalized;
 public:
  //! Constructor taking base name of file and
  //! and a header for the page which
  //! should be valid tex code. 
  EocOutput(std::string name,
            std::string descript) :
    outputFile(name)
  , level(0)
  , prevError(0)
  , prevSize(0)
  , initial(true)
  , finalized(false) {
    {
      std::ostringstream filestream;
      filestream << outputFile << ".tex";
      std::ofstream ofs(filestream.str().c_str(), std::ios::out);
      std::ostringstream filestreamBody;
      filestreamBody << outputFile << "_body.tex";
      std::ofstream ofsBody(filestreamBody.str().c_str(), std::ios::out);

      ofs << "\\documentclass[12pt,english]{article}\n"
         	<< "\\usepackage[T1]{fontenc}\n"
	        << "\\usepackage[latin1]{inputenc}\n"
	        << "\\usepackage{setspace}\n"
	        << "\\onehalfspacing\n"
	        << "\\makeatletter\n"
	        << "\\providecommand{\\boldsymbol}[1]{\\mbox{\\boldmath $#1$}}\n"
	        << "\\providecommand{\\tabularnewline}{\\\\}\n"
	        << "\\usepackage{babel}\n"
	        << "\\makeatother\n"
	        << "\\begin{document}\n"
          << "\\begin{center}\\large\n"
          << descript 
          << "\n\\end{center}\n\n"
          << "\\include{"
          << filestreamBody.str() 
          << "}\n";
      ofs.close();	
      ofsBody.close();	
    }
    {
      std::ostringstream filestream;
      filestream << outputFile << ".gnu";
      std::ofstream ofs(filestream.str().c_str(), std::ios::out);
      ofs.close();
    }
  }
  //! Destrucor writes the end document section to file
  //! if not done explicitly in the printTexEnd method
  ~EocOutput() {
    if (!finalized) {
      std::ostringstream filestream;
      filestream << outputFile << ".tex";
      std::ofstream ofs(filestream.str().c_str(), std::ios::app);
      std::ostringstream filestreamBody;
      filestreamBody << outputFile << "_body.tex";
      std::ofstream ofsBody(filestreamBody.str().c_str(), std::ios::app);
      ofsBody  << "\\end{tabular}\\\\\n\n";
	    ofs << "\\end{document}\n" << std::endl;
      ofs.close();
      ofsBody.close();
    }
  }
  //! write end of document information 
  //! and add the total runtime of the simulation
  void printTexEnd(double totaltime) {
    std::ostringstream filestream;
    filestream << outputFile << ".tex";
    std::ofstream ofs(filestream.str().c_str(), std::ios::app);
    std::ostringstream filestreamBody;
    filestreamBody << outputFile << "_body.tex";
    std::ofstream ofsBody(filestreamBody.str().c_str(), std::ios::app);
    ofsBody  << "\\end{tabular}\\\\\n\n"
	       << "Total time: " << totaltime << "\n";
	  ofs << "\\end{document}\n" << std::endl;
    ofs.close();
    ofsBody.close();
    finalized = true;
  }
  //! header information on setting for this simulation
  //! at the moment a grid class and 
  //! two other clases are passed to this method
  //! which should implement a method
  //!   printTexInfo(std::ostream&) const
  template <class GridType,class ProblemType,
            class SolverType>
  void printInput(const ProblemType& problem, 
                  const GridType& grid,
                  const SolverType& solver) {
    std::ostringstream filestreamBody;
    filestreamBody << outputFile << "_body.tex";
    std::ofstream ofsBody(filestreamBody.str().c_str(), std::ios::app);
    ofsBody  << "Grid: "      << grid.name() << "\\\\"
	           << "\n";
	  solver.printTexInfo(ofsBody);
    problem.printTexInfo(ofsBody);
    ofsBody.close();
  }
  //! \brief Write tabular for error and eoc output
  //!
  //! The eoc is computed using the size parameter
  //! and the error - previous size and errors are
  //! saved. The formular used is
  //! \f{eqnarray*} 
  //! \frac{\log(e_{\rm old}/e_{\rm new})}
  //!         {\log(S_{\rm old}/S_{\rm new})} 
  //! \f}
  //! where \f$ S \f$ is the size parameter. 
  //! 
  //! If for example the number of element \f$ N \f$
  //! in the grid is to be used to compute the eoc, then 
  //! \f{eqnarray*}
  //! S=\frac{1}{N^{\frac{1}{d}}}
  //! \f}
  //! where \f$ d \f$ is the dimension of the grid.
  //!
  //! \param size: the size of the grid which is used to compute the EOC
  //! \param error: a vector of computed error
  //! \param errDescr: a vector of column headers for the errors
  //! \param time: run time for this step of the simulation
  //! \param counter: could be time step counter or number of refinment
  //!                 steps for example
  void printTexAddError(double size,
                        std::vector<double> error, 
                        std::vector<std::string> errDescr,
                        double time, 
                        int counter) {
    std::ostringstream filestreamBody;
    filestreamBody << outputFile << "_body.tex";
    std::ofstream ofsBody(filestreamBody.str().c_str(), std::ios::app);
    std::ostringstream filestreamGnu;
    filestreamGnu << outputFile << ".gnu";
    std::ofstream ofsGnu(filestreamGnu.str().c_str(), std::ios::app);
		
    if (initial) {
      assert(errDescr.size()==error.size());
	    ofsBody << "\\begin{tabular}{|c|c|c|c|";
      for (int i=0;i<error.size();i++) {
        ofsBody << "|cc|";
      }
      ofsBody << "}\n"
	        << "\\hline \n"
          << "level & size & CPU-time & counter";
      for (int i=0;i<error.size();i++) {
        ofsBody << " & " << errDescr[i]
                << " & eoc ";
      }
	    ofsBody << "\n \\tabularnewline\n"
	            << "\\hline\n"
	            << "\\hline\n";
    } else {
      assert(error.size()==prevError.size());
    }
	  ofsBody <<  "\\hline \n"
	          << level << " & "
            << size << " & "
            << time << " & " 
            << counter;
	  ofsGnu << level << "  "
           << size << "  " 
           << time << "  " 
           << counter;
    for (int i=0;i<error.size();++i) {
      ofsBody << " & " << error[i] << " & ";
      ofsGnu << "    " << error[i] << " ";
      if (initial) {
        ofsBody << " --- ";
        ofsGnu << " 1/0 ";
      }
      else {
        double factor = prevSize/size;
        ofsBody << log(prevError[i]/error[i])/log(factor);
        ofsGnu << log(prevError[i]/error[i])/log(factor);
      }
    }
	  ofsBody << "\n"
            << "\\tabularnewline\n"
	          << "\\hline \n";
    ofsGnu << "\n";
    prevError = error;
    prevSize = size;
    level++;
    initial = false;
    ofsBody.close();
    ofsGnu.close();
  }

};

}
#endif
