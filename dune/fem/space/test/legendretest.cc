#include <config.h>

#include <iostream>

// #include <dune/grid/io/file/dgfparser/gridtype.hh>
// static const int dimw = dimworld;
#define S_GRID 0
#define YGRID 1
static const int dimw =2;
#if S_GRID 
#include<dune/grid/sgrid.hh> 

typedef Dune::SGrid<dimw,dimw> HGridType;
#endif

#if YGRID
#include<dune/grid/yaspgrid.hh> 

typedef Dune::YaspGrid< dimw > HGridType;
#endif

#include <dune/fem/operator/discreteoperatorimp.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/quadrature/cachingquadrature.hh> 

#include <dune/fem/gridpart/hierarchicgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh> 

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif
using namespace Dune;
using namespace Fem;

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
const int polOrd =1;/*POLORDER;*/

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

  //typedef Dune::SGrid<dimw,dimw> HGridType;
//! the index set we are using 
typedef LeafGridPart<HGridType> GridPartType;

//! define the function space, \f[ \R^2 \rightarrow \R \f]
// see dune/common/functionspace.hh
//typedef MatrixFunctionSpace < double , double, dimw , 3,5 > FuncSpace;
typedef FunctionSpace < double , double, dimw , 1 > FuncSpace;

//! define the function space our unkown belong to 
typedef  LegendreDiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
	polOrd,CachingStorage> DiscreteFunctionSpaceType;
//  typedef  DiscontinuousGalerkinSpace<FuncSpace, GridPartType, 
//  	polOrd,CachingStorage> DiscreteFunctionSpaceType;
//! define the type of discrete function we are using , see
typedef AdaptiveDiscreteFunction < DiscreteFunctionSpaceType > DiscreteFunctionType;

//! Get the Dofmanager type
typedef DofManager<HGridType> DofManagerType;


// the exact solution to the problem for EOC calculation 
struct ExactSolution
: public Fem::Function < FuncSpace , ExactSolution > 
{
  typedef FuncSpace::RangeType RangeType;
  typedef FuncSpace::JacobianRangeType JacobianRangeType;
  typedef FuncSpace::RangeFieldType RangeFieldType;
  typedef FuncSpace::DomainType DomainType;

  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate (const DomainType & x , RangeType & ret)  const
  {
    ret = 1.; // maximum of function is 2
     
       for(int i=0; i<DomainType::dimension; i++)
	ret *=sin(x[i]) ;
     }
  
  void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
  {
    evaluate ( x , ret );
  }

  void jacobian(const DomainType & x,JacobianRangeType & ret) const
  {
    
    for(int i=0; i<DomainType::dimension; i++)
      {
	double prod=1.0;
	for(int j=0;j<DomainType::dimension;j++){
	  if(j==i)
	    prod*=cos(x[j]);
	  else
	    prod*=sin(x[j]);
	}
       
	ret[0][i]=prod;
      }
    
  }




};
 
// ********************************************************************
template <class DiscreteFunctionType>
class L2Projection
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

 public:
  template <class FunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc, int polOrd) 
  {
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType Iterator;

    const DiscreteFunctionSpaceType& space =  discFunc.space();

    discFunc.clear();

    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    RangeType ret (0.0);
    std::vector< RangeType > phi;
    //diagomal of massmatrix
    DiscreteFunctionType mass("mass",space);
    mass.clear();

    Iterator endit = space.end();
    for(Iterator it = space.begin(); it != endit ; ++it) 
    {
      // Get quadrature rule
      CachingQuadrature<GridPartType,0> quad(*it, polOrd);

      LocalFuncType lf = discFunc.localFunction(*it);
      LocalFuncType tmp = mass.localFunction(*it);
      
      //! Note: basis functions must be ortho-normal!!!!
      typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType ; 
      const BasisFunctionSetType & basisSet = lf.basisFunctionSet();

      const typename HGridType::template Codim<0>::Entity::Geometry& 
        itGeom = (*it).geometry();
      
      const int quadNop = quad.nop();
      const int numDofs = lf.numDofs();
      phi.resize( numDofs );

      for(int qP = 0; qP < quadNop ; ++qP) 
      {// double det = (*it).geometry().integrationElement( quad.point(qP) );
        f.evaluate(itGeom.global(quad.point(qP)), ret);

        basisSet.evaluateAll( quad[ qP ], phi );

        for(int i=0; i<numDofs; ++i) 
        {
	  //	  tmp[i]+=quad.weight(qP)*SQR(phi)*det ;
          lf[i] += quad.weight(qP) * (ret * phi[i])/*det*/ ;
        }
      }
//       for(int i=0; i<numDofs; ++i) {
// 	lf[i]/=tmp[i];}
   }
  }
  
  template <class FunctionType>
  static void project (const FunctionType &f, DiscreteFunctionType &discFunc) 
  {
    const DiscreteFunctionSpaceType& space =  discFunc.space();
    int polOrd = 2 * space.order();
    project(f,discFunc,polOrd);
  }
};


// calculates || u-u_h ||_H1
template <class DiscreteFunctionType>
class L2Error
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: RangeType RangeType;
  typedef typename DiscreteFunctionSpaceType :: JacobianRangeType JacobianRangeType;

public:
  template <class FunctionType>
  RangeType h1norm (const FunctionType &f, DiscreteFunctionType &discFunc,
      double time, int polOrd) const
  {
    const DiscreteFunctionSpaceType & space = discFunc.space();

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    RangeType ret (0.0);
    RangeType phi (0.0);
    RangeType tmp (0.0);

    JacobianRangeType psi(0.0);
    JacobianRangeType xi(0.0);
 


    RangeType error(0.0);

    enum { dimRange = DiscreteFunctionSpaceType :: DimRange };
    enum { dimDomain =DiscreteFunctionSpaceType :: DimDomain};
    IteratorType endit = space.end();
    for(IteratorType it = space.begin(); it != endit ; ++it)
    {
      CachingQuadrature<GridPartType,0> quad(*it, polOrd);
      LocalFuncType lf = discFunc.localFunction(*it);
      const int quadNop = quad.nop();
      for(int qP = 0; qP < quadNop; ++qP)
      {
        double weight = quad.weight(qP) * (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
	f.jacobian((*it).geometry().global(quad.point(qP)),psi);
        lf.evaluate(quad[qP],phi);
	lf.jacobian(quad[qP],xi);
	
	tmp=0.0;

	for(int i=0; i< dimDomain; ++i)
	  {
	    for(int j=0;j< dimRange; ++j)
	      
	      tmp[j]+=SQR(xi[j][i]-psi[j][i]);	
	  }

	for(int i=0; i< dimRange; ++i)
	  error[i] += weight * (SQR(ret[i] - phi[i])+tmp[i]);
      }
    }
    
    
    for(int i=0; i< dimRange; ++i) 
      error[i] = sqrt(error[i]);
    
    return error;
  }

template <class FunctionType>
  RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
      double time, int polOrd) const
  {
    const DiscreteFunctionSpaceType & space = discFunc.space();

    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    typedef typename DiscreteFunctionType::LocalFunctionType LocalFuncType;

    RangeType ret (0.0);
    RangeType phi (0.0);

    RangeType error(0.0);

    enum { dimRange = DiscreteFunctionSpaceType :: DimRange };

    IteratorType endit = space.end();
    for(IteratorType it = space.begin(); it != endit ; ++it)
    {
      CachingQuadrature<GridPartType,0> quad(*it, polOrd);
      LocalFuncType lf = discFunc.localFunction(*it);
      const int quadNop = quad.nop();
      for(int qP = 0; qP < quadNop; ++qP)
      {
        double weight = quad.weight(qP) * (*it).geometry().integrationElement(quad.point(qP));
        f.evaluate((*it).geometry().global(quad.point(qP)),time, ret);
       
	lf.evaluate(quad[qP],phi);

        for(int i=0; i< dimRange; ++i)
          error[i] += weight * SQR(ret[i] - phi[i]);
      }
    }
    
    for(int i=0; i< dimRange; ++i) 
      error[i] = sqrt(error[i]);
    
    return error;
  }



  template <class FunctionType>
  RangeType h1norm (const FunctionType &f, DiscreteFunctionType &discFunc,
		    double time) const
  {
    const DiscreteFunctionSpaceType & space = discFunc.space();
    int polOrd = 2 * space.order() + 2;
    return h1norm(f,discFunc,time,polOrd);
  }









  template <class FunctionType>
  RangeType norm (const FunctionType &f, DiscreteFunctionType &discFunc,
      double time) const
  {
    const DiscreteFunctionSpaceType & space = discFunc.space();
    int polOrd = 2 * space.order() + 2;
    return norm(f,discFunc,time,polOrd);
  }
};
// ********************************************************************
double algorithm (HGridType& grid, DiscreteFunctionType& solution  , int turn )
{
   GridPartType part ( grid );
   DiscreteFunctionSpaceType linFuncSpace ( part );
   ExactSolution f;
   L2Error < DiscreteFunctionType > l2err;
       
   //! perform l2-projection
   L2Projection<DiscreteFunctionType>::
     project(f, solution);

   // calculation L2 error 
   // pol ord for calculation the error chould by higher than 
   // pol for evaluation the basefunctions 
   typedef DiscreteFunctionSpaceType :: RangeType RangeType; 
   RangeType error = l2err.norm(f ,solution, 0.0);
   RangeType h1error= l2err.h1norm(f,solution,0.0);
   for(int i=0; i<RangeType::dimension; ++i)
     {
       std::cout << "\nL2 Error["<<i<<"] : " << error[i] << "\n\n";
        std::cout << "\nH1 Error["<<i<<"] : " << h1error[i] << "\n\n";
     }
#if HAVE_GRAPE
   // if Grape was found, then display last solution 
   if(0 && turn > 0)
   {
     GrapeDataDisplay < HGridType > grape(part); 
     grape.dataDisplay( solution );
   }
#endif
   
   //return sqrt(error*error);
   return sqrt(h1error*h1error);
}


//**************************************************
//
//  main programm, run algorithm twice to calc EOC 
//
//**************************************************
int main (int argc, char **argv)
{
  MPIManager::initialize( argc, argv );

  if(argc != 2)
  {
    fprintf(stderr,"usage: %s <maxlevel> \n",argv[0]);
    exit(1);
  }

  int ml = atoi( argv[1] );
  double* error = new double[ml];
  // char tmp[16]; sprintf(tmp,"%d",dimw);
//   std::string macroGridName (tmp); 
//   macroGridName += "dgrid.dgf";

#if S_GRID
  Dune::FieldVector<int,dimw> N(2);
  Dune::FieldVector<HGridType::ctype,dimw> L(0.0);
  Dune::FieldVector<HGridType::ctype,dimw> H(1.0);
  HGridType grid(N,L,H);
#endif

#if YGRID
FieldVector<double, dimw> lang;
 for(int i=0;i<dimw;i++)
   lang[i]= 1.0;
  
  FieldVector<int, dimw> anz;
  for(int i=0;i<dimw;i++)
    anz[i] = 1;
  FieldVector<bool, dimw> per;
  for(int i=0;i<dimw;i++)
    per[i] = false; 
  
  HGridType grid(lang,anz,per,1);

#endif

  
 //  GridPtr<HGridType> gridptr(macroGridName);
//   HGridType& grid=*gridptr;
 //  const int step = Dune::DGFGridInfo<HGridType>::refineStepsForHalf();

  const int step =1;

  GridPartType part ( grid );
  DiscreteFunctionSpaceType linFuncSpace ( part );
  DiscreteFunctionType solution ( "sol", linFuncSpace );
  solution.clear();
  
  for(int i=0; i<ml; i+=step)
  {
    grid.globalRefine(step);
    DofManagerType& dm = DofManagerType :: instance( grid );
    dm.resize();
    error[i] = algorithm ( grid , solution , i==ml-1);
    if (i>0) {
      double eoc = log( error[i-step]/error[i]) / M_LN2; 
      std::cout << "EOC = " << eoc << " \n";
   
    }
  }
  delete [] error;
  
  return 0;
}

