#undef HAVE_MPI
#include "../../../macrogridparser/gridtype.hh"
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/albertagrid.hh>
#include "advectdiff.hh"
#include "odesolver.hh"


// Modelle
// #define POLORDER 1
// Approximations Ordnung
enum {order=POLORDER,rksteps=POLORDER+1}; 
// Gitterauswahl
// typedef YaspGrid<2,2> GridType;
// typedef SGrid<2,2> GridType;
// typedef AlbertaGrid<2,2> GridType;

// Modell- und Flussauswahl
// Skalar
#if PROBLEM = 1
   #include "scalarmodels.hh"
   typedef U0<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   InitialDataType problem(0.01,true);
   typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplTimeStepper<DgType> ODEType;
#slif PROBLEM = 2
// Euler
   #include "euler_mhd/eulermodel.hh"
   typedef U0RotatingCone InitialDataType;
   typedef EulerModel<GridType,InitialDataType> ModelType;
   typedef DWNumFlux<ModelType> FluxType;
   InitialDataType problem;
   typedef DGAdvectionOperator<ModelType,DWNumFlux,order> DgType;
   typedef DuneODE::ExplTimeStepper<DgType> ODEType;
#endif
// *** Operator typedefs
// Timestepper:
// typedef DuneODE::ExplRungeKutta<DgType> ODEType;
