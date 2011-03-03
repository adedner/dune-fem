#ifndef DUNE_FEM_IO_HH
#define DUNE_FEM_IO_HH

#include <iostream>

namespace Dune
{

  /** \brief create a directory
   *
   *  \param[in]  name  name of the directory to create
   *
   *  \returns whether the directory has been successfully created
   */
  bool createDirectory ( const std::string &name );

  /** \brief check whether a directory exists 
   *
   *  \param[in]  name  name of the directory to create
   *
   *  \returns true if directory exists, false otherwise 
   */
  bool directoryExists ( const std::string &name );
}

#endif // #ifndef DUNE_FEM_IO_HH
