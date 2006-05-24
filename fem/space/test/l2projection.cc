#ifndef DIM_OF_WORLD 
static const int dimw = 2;
#else
static const int dimw = DIM_OF_WORLD;
#endif

#ifndef DIM
static const int dimp = 2;
#else
static const int dimp = DIM;
#endif

#include <iostream>
#include <config.h>

#define SGRID 0
#define AGRID 1

#include <../../operator/discreteoperatorimp.hh>
#include <../lagrangespace.hh>
#include <../../discretefunction/dfadapt.hh>
#include "../dgspace.hh"
#include "../../quadrature/cachequad.hh"

#include "../leafindexset.hh"
#include <dune/grid/common/gridpart.hh>

#include <dune/grid/common/referenceelements.hh>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

using namespace Dune;

#if SGRID
#include <dune/grid/sgrid.hh>
typedef SGrid  < dimp, dimw > GridType;
#endif

#if AGRID  
#include <dune/grid/albertagrid.hh>
typedef AlbertaGrid< dimp, dimw > GridType;
#endif

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
const int polOrd = POLORDER;

//***********************************************************************
/*! L2 Projection of a function f: 

  This is an example how to solve the equation on 
  \f[\Omega = (0,1)^2 \f]

  \f[ \int_{\Omega} u \phi = \int_{\Omega} f \phi  \ \ \ in \Omega \f]
  \f[ f(x,y) = x ( 1 - x) y ( 1 - y ) \f]

  Here u is the L_2 projection of f. 

  The Projection should converge to the given function f.
  with the finite element method using lagrangian elements of polynom order +1.
*/
//***********************************************************************

//! the index set we are using 
typedef HierarchicGridPart<GridType> GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
typedef FunctionSpace < double , double, dimp , 1 > FuncSpace;

//! define the function space our unkown belong to 
//! see dune/fem/lagrangebase.hh
typedef DiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
  polOrd,CachingStorage> DiscreteFunctionSpaceType;

//! define the type of discrete function we are using , see
//! dune/fem/discfuncarray.hh
typedef DFAdapt < DiscreteFunctionSpaceType > DiscreteFunctionType;

//! the exact solution to the problem for EOC calculation 
class ExactSolution : public Function < FuncSpace , ExactSolution > 
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;
public:
  ExactSolution (FuncSpace &f) : Function < FuncSpace , ExactSolution > ( f ) {}
 
  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret = 2.; // maximum of function is 2
    for(int i=0; i<DomainType::dimension; i++)
      ret *= x[i]*(1.0 -x[i])*4.;
  }
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }
};
 
// ********************************************************************
template <class DiscreteFunctionType>
class L2Projection
{
  typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;

 public:
  template <class FunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc, int polOrd) 
  {
    typedef typename DiscreteFunctionSpaceType::Traits::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::Traits::IteratorType Iterator;

    const DiscreteFunctionSpaceType& space =  discFunc.getFunctionSpace();

    discFunc.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typename DiscreteFunctionSpaceType::RangeType ret (0.0);
    typename DiscreteFunctionSpaceType::RangeType phi (0.0);

    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      // Get quadrature rule
      CachingQuadrature<GridType,0> quad(*it, polOrd);

      LocalFuncType lf = discFunc.localFunction(*it);

      //! Note: BaseFunctions must be ortho-normal!!!!
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType ; 
      const BaseFunctionSetType & baseset =
        lf.getBaseFunctionSet();

      const typename GridType::template Codim<0>::Entity::Geometry& 
        itGeom = (*it).geometry();
     
      const int quadNop = quad.nop();
      const int numDofs = lf.numDofs();
      for(int qP = 0; qP < quadNop ; ++qP) 
      {
        f.evaluate(itGeom.global(quad.point(qP)), ret);
        for(int i=0; i<numDofs; ++i) {
          baseset.eval(i,quad,qP,phi);
          lf[i] += quad.weight(qP) * (ret * phi) ;
        }
      }
    }
  }
  
  template <class FunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) 
  {
    const DiscreteFunctionSpaceType& space =  discFunc.getFunctionSpace();
    int polOrd = 2 * space.polynomOrder();
    project(f,discFunc,polOrd);
  }
};


// calculates || u-u_h ||_L2
template <class DiscreteFunctionType>
class L2Error
{
  typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;

public:
  template <class FunctionType>
  RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
      double time, int polOrd) const
  {
    const DiscreteFunctionSpaceType & space = discFunc.getFunctionSpace();

    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    RangeType ret (0.0);
    RangeType phi (0.0);

    RangeType error(0.0);

    enum { dimRange = DiscreteFunctionSpaceType :: DimRange };

    IteratorType endit = space.end();
    for(IteratorType it = space.begin(); it != endit ; ++it)
    {
      CachingQuadrature<GridType,0> quad(*it, polOrd);
      LocalFuncType lf = discFunc.localFunction(*it);
      const int quadNop = quad.nop();
      for(int qP = 0; qP < quadNop; ++qP)
      {
        double weight = quad.weight(qP) * (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
        lf.evaluate((*it),quad,qP,phi);

        for(int i=0; i< dimRange; ++i)
          error[i] += weight * SQR(ret[i] - phi[i]);
      }
    }

    for(int i=0; i< dimRange; ++i) 
      error[i] = sqrt(error[i]);
    
    return error;
  }

  template <class FunctionType>
  RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
      double time) const
  {
    const DiscreteFunctionSpaceType & space = discFunc.getFunctionSpace();
    int polOrd = 2 * space.polynomOrder() + 2;
    return norm(f,discFunc,time,polOrd);
  }
};
// ********************************************************************
double algorithm (GridType& grid, int turn )
{
   GridPartType part ( grid );

   DiscreteFunctionSpaceType linFuncSpace ( part );
   DiscreteFunctionType solution ( "sol", linFuncSpace );
   solution.clear();
      
   ExactSolution f ( linFuncSpace ); 
   L2Error < DiscreteFunctionType > l2err;
       
   //! perform l2-projection
   L2Projection<DiscreteFunctionType>::
     project(f, solution);

   // calculation L2 error 
   // pol ord for calculation the error chould by higher than 
   // pol for evaluation the basefunctions 
   typedef DiscreteFunctionSpaceType :: RangeType RangeType; 
   RangeType error = l2err.norm(f ,solution, 0.0);
   for(int i=0; i<RangeType::dimension; ++i)
     std::cout << "\nL2 Error["<<i<<"] : " << error[i] << "\n\n";
  
#if HAVE_GRAPE
   // if Grape was found, then display last solution 
   if(turn > 0)
   {
     GrapeDataDisplay < GridType > grape(grid); 
     grape.dataDisplay( solution );
   }
#endif
   
   return error;
}


//**************************************************
//
//  main programm, run algorithm twice to calc EOC 
//
//**************************************************
int main (int argc, char **argv)
{
  if(argc != 2)
  {
    fprintf(stderr,"usage: %s <maxlevel> \n",argv[0]);
    exit(1);
  }
  int ml = atoi( argv[1] );
  double* error = new double[ml];
  char tmp[16]; sprintf(tmp,"%d",dimp);
  std::string macroGridName (tmp); 
  macroGridName += "dgrid.al";

#if SGRID 
  const int step = 1;
#else 
  const int step = 2;
#endif
#if SGRID
   // this leads to the same number of points for SGrid and AlbertGrid
   int n[dimp];
   double h[dimp];
   for(int i=0; i<dimp; i++)  { n[i] = 2; h[i] = 1.0; }

   GridType grid ((int *) &n, (double *) &h );
#else
   GridType grid ( macroGridName.c_str() );
#endif

  
  for(int i=0; i<ml; i+=step)
  {
    grid.globalRefine(step);
    error[i] = algorithm ( grid , 0);
    if (i>0) {
      double eoc = log( error[i-step]/error[i]) / M_LN2; 
      std::cout << "EOC = " << eoc << " \n";
    }
  }
  delete [] error;
  return 0;
}

