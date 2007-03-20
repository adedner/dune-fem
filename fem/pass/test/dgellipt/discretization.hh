/***********************************************************************************************
 Discretization is the interface class used to define the discretization of the problem
***********************************************************************************************/
#ifndef __LDGEXAMPLE_DISCRETIZATION_HH__
#define __LDGEXAMPLE_DISCRETIZATION_HH__


/* include definition of the physical problem (coefficient functions, etc.) */
#include "model.hh"

// choose discrete function type
#include <dune/fem/space/lagrangespace.hh>

#include <dune/grid/common/gridpart.hh>

// ascii parser for parameter files 
#include <dune/fem/io/file/asciiparser.hh>

// grape data io 
#include <dune/fem/io/file/grapedataio.hh>

// if grape was configured then include headers 
#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/grid/io/visual/combinedgrapedisplay.hh>
#endif

#include "discretemodels.hh"

#include <dune/fem/operator/inverseoperators.hh>
#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>

#include <dune/fem/pass/dgelliptpass.hh>

// definition of L2Error 
#include <dune/fem/misc/l2error.hh>

// implementation of ldg flux 
#include <dune/fem/pass/ldgflux.hh>

#include <dune/fem/space/dgspace/dgadaptmanager.hh>

using namespace Dune;

namespace LDGExample { 

template <class ModelImpType, int polOrd=0 >
struct DiscrParam
{
  typedef ModelImpType  ModelType;

  enum { dimRange = ModelType::dimRange };
 
  typedef typename ModelType::GridType             GridType;

  enum { dim      = GridType::dimension };
  enum { dimworld = GridType::dimensionworld };
  enum { polyOrder = polOrd };

  typedef typename ModelType::FieldType            FieldType;
  typedef typename ModelType::FuncSpaceType        FuncSpaceType;

  typedef LeafGridPart < GridType > GridPartType ;
};


// The actual operator
template <class GradientModelType, class LaplaceModelType,
          class VelocityModelType> 
class MySpaceOperator :
 public Operator<
    double,
    double,
    typename VelocityModelType::Traits::DestinationType,
    typename VelocityModelType::Traits::DestinationType>
{
  typedef typename VelocityModelType :: Traits Traits;
public:
  typedef MySpaceOperator<GradientModelType,
    LaplaceModelType,VelocityModelType> ThisType;
  enum { polOrd = LaplaceModelType :: polynomialOrder };
  
  typedef typename Traits:: DestinationType DestinationType;
  typedef typename Traits:: DestinationType DiscreteFunctionType;

  typedef typename Traits::GridPartType GridPartType ;
  typedef typename GridPartType :: Traits :: GridType GridType;

  typedef typename LaplaceModelType::Traits::DestinationType SolutionType;
  typedef typename SolutionType :: DiscreteFunctionSpaceType LastSpaceType;

  typedef DofManager<GridType> DofManagerType; 
  typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

  typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename GradientModelType::Traits::DiscreteFunctionSpaceType GradDiscreteFunctionSpaceType;

  typedef DiscreteFunctionSpaceType VeloSpaceType;
  typedef StartPass<DiscreteFunctionType> Pass0Type;
  // note, the destination type of the pass 0 is the argument type of pass 1
  typedef LocalDGElliptGradientPass<GradientModelType , Pass0Type> GradPassType;
  typedef LocalDGElliptPass<LaplaceModelType, GradPassType> LastPassType;
  typedef LocalDGPass<VelocityModelType,LastPassType> VeloPassType;

  typedef GradDiscreteFunctionSpaceType GradSpaceType;
  
  typedef typename LastSpaceType :: FunctionSpaceType FuncSpaceType;
  typedef typename FuncSpaceType::RangeType RangeType;

  //! the exact solution to the problem for EOC calculation 
  class ExactSolution : public Function < FuncSpaceType , ExactSolution >
  {
    typedef typename FuncSpaceType::RangeType RangeType;
    typedef typename FuncSpaceType::RangeFieldType RangeFieldType;
    typedef typename FuncSpaceType::DomainType DomainType;
  public:
    ExactSolution (FuncSpaceType &f) : Function < FuncSpaceType , ExactSolution > ( f ) {}

    //! see problem.cc 
    void evaluate (const DomainType & x , RangeType & ret) const
    {
      ret = exactSolution( &x[0] );
    }
    void evaluate (const DomainType & x , RangeFieldType time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };

  //! the exact solution to the problem for EOC calculation 
  template <class FuncSPCType>
  class ExactGradient : public Function < FuncSPCType , ExactGradient<
                        FuncSPCType > >
  {
    typedef typename FuncSPCType::RangeType RangeType;
    typedef typename FuncSPCType::RangeFieldType RangeFieldType;
    typedef typename FuncSPCType::DomainType DomainType;
  public:
    ExactGradient (const FuncSPCType &f) 
      : Function < FuncSPCType , ExactGradient< FuncSPCType > > ( f ) {}

    void evaluate (const DomainType & x , RangeType & ret) const
    {
      exactGradient( &x[0], &ret[0] );
    }

    void evaluate (const DomainType & x , double time , RangeType & ret) const
    {
      evaluate ( x , ret );
    }
  };

  MySpaceOperator ( GridType & grid , 
                    GradientModelType & gm, 
                    LaplaceModelType & lpm , 
                    VelocityModelType& vm,
                    std::string paramfile)
    : grid_(grid)
    , dm_(DofManagerFactoryType::getDofManager(grid_))
    , gridPart_(grid_)
    , gradSpace_(gridPart_)
    , lastSpace_(gridPart_)
    , veloSpace_(gridPart_)
    , pass0_()
    , pass1_( gm , pass0_, gradSpace_)
    , lastPass_( lpm, pass1_, lastSpace_ , paramfile )
    , veloPass_( vm, lastPass_, veloSpace_ )
    , steps_(2)
  {
    readParameter(paramfile,"EOCSteps",steps_);

    std::cout << "Created SpaceOperator \n";
    std::cout.flush();
  }

  void testConsecutive() 
  {
    typedef typename LastSpaceType :: IndexSetType IndexSetType; 
    const IndexSetType& iset = lastSpace_.indexSet();
    
    typedef typename LastSpaceType :: IteratorType IteratorType;
    std::vector<bool> visited(iset.size(0),false);
    IteratorType endit = lastSpace_.end();
    size_t count = 0;
    for(IteratorType it = lastSpace_.begin(); it != endit ; ++it) 
    {
      ++count;
      int idx = iset.index(*it);
      visited[idx] = true;
    }

    assert( count == visited.size() );
    for(size_t i=0; i<visited.size(); ++i) 
    {
      assert( visited[i] == true );
    }
  }

  void operator()(const DestinationType& arg, DestinationType& dest) const 
  {
    const_cast<ThisType&> (*this).apply(arg,dest);
  }

  void adaptGrid (DestinationType& dest) 
  {
    typedef typename LastPassType :: RestrictProlongOperatorType
      RPOpType;

    typedef AdaptationManager <GridType, RPOpType> AdaptManagerType;
    AdaptManagerType adop(grid_,lastPass_.restrictProlongOperator() );

    typedef typename LastSpaceType :: IteratorType IteratorType;
    std::cout << "Old size of space is " << lastSpace_.size() << "\n";
    int count = 0;
    
    typedef typename DestinationType :: DomainType DomainType; 
    DomainType point(0.5);
  
    //int halfSize = lastSpace_.indexSet().size(0)/2;
    IteratorType endit = lastSpace_.end();
    for(IteratorType it = lastSpace_.begin(); it != endit ; ++it) 
    {
      if( (it->geometry()[0] - point).two_norm() < 0.3) 
      //if(count < halfSize) 
      //if(count % 2 == 0) 
      {
        //std::cout << "Mark entity \n";
        grid_.mark(3,it);
      }
      ++count;
    }

    adop.adapt();
    std::cout << "New size of space is " << lastSpace_.size() << "\n";

    testConsecutive();
  }
  
  // apply space discretisation 
  void apply(const DestinationType& arg, DestinationType& velo)
  {
    typedef typename GradientModelType::Traits::DestinationType  GradFuncType; 
    typedef typename GradientModelType::Traits::DiscreteFunctionSpaceType GradientDiscreteFunctionSpaceType; 
    typedef typename GradientDiscreteFunctionSpaceType :: RangeType  GradRangeType;

    std::vector<RangeType> error(steps_);
    std::vector<GradRangeType> gradError(steps_);
    
    //DestinationType & Arg = const_cast<DestinationType&> (arg);
    //adaptGrid(Arg);

    DestinationType & Arg = const_cast<DestinationType&> (arg);
    adaptGrid(Arg);

    for(int i=0; i<steps_; ++i)
    {
      if(i > 0)
      {
        // refineGlobal is defined in description.hh
        grid_.globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());
        dm_.resize();
        dm_.dofCompress();
      }

      FuncSpaceType sp; 
      ExactSolution exact(sp); 

      SolutionType& dest = const_cast<SolutionType&> (lastPass_.destination());
      dest.clear();

      veloPass_(arg,velo);
      {
        // only for LDG method
        //velo.clear();
        //lastPass_.evalGradient(dest,velo);
        
        L2Error < DestinationType > l2errGrad;
        ExactGradient< VeloSpaceType > exactGrad(gradSpace_);

        gradError[i] = l2errGrad.norm(exactGrad , velo);
      }
      
#if HAVE_GRAPE
      GrapeDataDisplay < GridType > grape( gridPart_.grid() ); 
      grape.addData( velo );
      grape.dataDisplay( dest );
#endif

      L2Error < SolutionType > l2err;
      // pol ord for calculation the error chould by higher than 
      // pol for evaluation the basefunctions 
      error[i] = l2err.norm(exact , dest);

      for(int k=0; k<RangeType::dimension; k++)
        std::cout << "\nError["<<k<<"] : " << error[i][k] << "\n";
      
      for(int k=0; k<GradRangeType::dimension; k++)
        std::cout << "GradError["<<k<<"] : " << gradError[i][k] << "\n";
    
      if( i > 0 )
      {
        for(int k=0; k<RangeType::dimension; k++)
        {
          double eoc = log( error[i-1][k]/error[i][k]) / M_LN2;
          std::cout << "\nEOC["<<i <<"] = " << eoc << " \n";
        }
        
        for(int k=0; k<GradRangeType::dimension; k++)
        {
          double eoc = log( gradError[i-1][k]/gradError[i][k]) / M_LN2;
          std::cout << "Grad EOC["<<i <<"] = " << eoc << " \n";
        }
      }
    }

  }

  SolutionType& solution () 
  {
    SolutionType& dest = const_cast<SolutionType&> (lastPass_.destination());
    return dest;
  }
  
  template <class TimeProviderType>
  void timeProvider(TimeProviderType & tp )
  {
    pass1_.timeProvider(&tp);
  }

  GridPartType & gridPart () { return gridPart_; }

  DestinationType * createDestinationFct (std::string name) 
  {
    return new DestinationType (name , veloSpace_ );
  }

private:
  GridType & grid_;
  DofManagerType & dm_;

  // we use the same index set and grid part for all spaces
  GridPartType gridPart_;

  mutable GradSpaceType gradSpace_;
  mutable LastSpaceType lastSpace_;
  mutable VeloSpaceType veloSpace_;

  mutable Pass0Type pass0_;
  mutable GradPassType pass1_;
  mutable LastPassType lastPass_;
  mutable VeloPassType veloPass_;

  int steps_; 
};

template <class DiscrType> 
void simul(typename DiscrType::ModelType & model, std::string paramFile) 
{
  typedef typename DiscrType::GridType                     GridType;
  typedef typename DiscrType::ModelType             ModelType;
  enum { polOrd = DiscrType::polyOrder };

  // choice of fluxes 
  typedef LDGFlux<ModelType> NumericalFluxType;

  typedef LaplaceDiscreteModel < ModelType, NumericalFluxType, polOrd > LaplaceModelType;
  typedef GradientDiscreteModel < ModelType, NumericalFluxType, polOrd-1 > GradientModelType;
  typedef VelocityDiscreteModel < ModelType, polOrd-1 > VelocityModelType;
  
  typedef MySpaceOperator <  GradientModelType, 
                             LaplaceModelType,
                             VelocityModelType> 
                SpaceOperatorType; 
  typedef typename SpaceOperatorType :: DestinationType DestinationType;
   
  // initialize grid
  char dummyfile [4096];
  const char * paramfile = paramFile.c_str();
  readParameter(paramfile,"Grid",dummyfile);
  std::string macroGridName(dummyfile);

  int level;
  readParameter(paramfile,"StartLevel",level);

  GridPtr<GridType> gridptr(macroGridName,MPIHelper::getCommunicator()); 
  GridType & grid = *gridptr;

  grid.globalRefine(DGFGridInfo<GridType>::refineStepsForHalf() * level);
  grid.loadBalance();

  //int bla;
  //fscanf(stdin,"%d",&bla);
  
  std::cout << "Grid size = " << grid.size(0) << "\n";

  int display = 0;
  readParameter(paramfile,"display",display);

  // read parameter for LDGFlux 
  double beta = 0.0;
  if(!readParameter(paramfile,"beta",beta))
  {
    std::cout << "Using beta = "<< beta << "\n";
  }
  double power = 1.0;
  if(!readParameter(paramfile,"power",power))
  {
    std::cout << "Using power of h = "<< power << "\n";
  }
  double eta = 1.0;
  if(!readParameter(paramfile,"eta",eta))
  {
    std::cout << "Using eta = "<< eta << "\n";
  }
  NumericalFluxType numericalFlux(model,beta,power,eta);
  
  LaplaceModelType lpm(model, numericalFlux);
  GradientModelType gm(model, numericalFlux);
  VelocityModelType vm(model);

  SpaceOperatorType spaceOp(grid , gm, lpm , vm, paramfile );
  //fscanf(stdin,"%d",&bla);
  
  //! storage for the discrete solution and its update
  DestinationType *solution = spaceOp.createDestinationFct("solution");
  DestinationType *tmpRhs = spaceOp.createDestinationFct("tmpRhs");

  // initial data != 0
  solution->clear();
  spaceOp(*tmpRhs,*solution);

#if HAVE_GRAPE
  if( display == 1 )
  {
    GrapeDataDisplay < GridType > grape( spaceOp.gridPart() ); 
    grape.dataDisplay( spaceOp.solution() );
  }
#endif
  delete solution; 
  delete tmpRhs;
}

} // end namespace LDGExample
#endif
