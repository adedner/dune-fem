//- system includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstring>
#include <pthread.h>
#include <cassert>
#include <sys/times.h>
#include <utility>

#include <config.h>

namespace parDG_NoMPI {

#include "emptycommunicator.hpp"
#define COMMUNICATOR_HPP

static const int MPI_SUM = 0;

#include "bulirsch_stoer.cpp"
#include "iterative_solver.cpp"
#include "ode_solver.cpp"
#include "sirk.cpp"

#include "matrix.cpp"
#include "qr_solver.cpp"
#include "lu_solver.cpp"
#include "ssp.cpp"
#include "dirk.cpp"
#include "runge_kutta.cpp"
#include "gmres.cpp"
#include "fgmres.cpp"
#include "bicgstab.cpp"
#include "cg.cpp"
#include "vector.cpp"

#include "quadrature0d.cpp"
#include "quadrature1d.cpp"
#include "quadrature2d.cpp"
#include "quadrature3d.cpp"

} // end namespace parDG_MPI
