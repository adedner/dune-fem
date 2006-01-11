#ifndef DUNE_DGPASS_HH
#define DUNE_DGPASS_HH

#include "pass.hh"
#include "selection.hh"
#include "discretemodel.hh"
#include "modelcaller.hh"

// * needs to move
// #include "../misc/timenew.hh"
#include "../misc/timeutility.hh"

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/utility/twistutility.hh>


namespace Dune {

  //! Concrete implementation of Pass for LDG.
  template <class DiscreteModelImp, class PreviousPassImp>
  class LocalDGPass :
    public LocalPass<DiscreteModelImp, PreviousPassImp> 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    // Types extracted from the underlying grids
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;
    
    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    LocalDGPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                DiscreteFunctionSpaceType& spc,
		int quadOrd0=-1,int quadOrd1=-1) :
      BaseType(pass, spc),
      caller_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      dtMin_(std::numeric_limits<double>::max()),
      fMat_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0),
      source_(0.0),
      grads_(0.0),
      time_(0),
      diffVar_(),
      twistUtil_(spc.grid()),
      quadOrd0_(quadOrd0),
      quadOrd1_(quadOrd1)
    {
      if (quadOrd0_==-1)
	quadOrd0_ = 2*spc_.polynomOrder();
      if (quadOrd1_==-1)
	quadOrd1_ = 2*spc_.polynomOrder()+1;
    }
   
    //! Destructor
    virtual ~LocalDGPass() {
      //delete caller_;
    }

    //! Stores the time provider passed by the base class in order to have
    //! access to the global time
    virtual void processTimeProvider(TimeProvider* time) {
      time_ = time;
    }

    //! Estimate for the timestep size
    double timeStepEstimate() const {
      return dtMin_;
    }

  private:
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      dest_->clear();

      caller_.setArgument(*arg_);

      // time initialisation
      dtMin_ = std::numeric_limits<double>::max();
      if (time_) {
        caller_.setTime(time_->time());
      }
      else {
        caller_.setTime(0.0);
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      if (time_) {
        time_->provideTimeStepEstimate(dtMin_);
      }

      caller_.finalize();
    }

    void applyLocal(EntityType& en) const
    {
      //- typedefs
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
      

      //- statements
      caller_.setEntity(en);
      LocalFunctionType updEn = dest_->localFunction(en);
      //GeometryType geom = en.geometry().type();
      
      VolumeQuadratureType volQuad(en, quadOrd0_);

      double vol = volumeElement(en, volQuad);
      //std::cout << "Vol = " << vol << std::endl;
      
      const IndexSetType& iset = spc_.indexSet();
      const BaseFunctionSetType& bsetEn = spc_.getBaseFunctionSet(en);

      // Volumetric integral part
      for (int l = 0; l < volQuad.nop(); ++l) {
        caller_.analyticalFlux(en, volQuad, l, fMat_);
        caller_.source(en, volQuad, l, source_);

        for (int i = 0; i < updEn.numDofs(); ++i) {
          updEn[i] += 
            (bsetEn.evaluateGradientSingle(i, en, volQuad, l, fMat_) +
             bsetEn.evaluateSingle(i, volQuad, l, source_))*
            volQuad.weight(l)*
            en.geometry().integrationElement(volQuad.point(l))/vol;
        }
      }
     
      // Surface integral part
      IntersectionIterator endnit = en.iend();
      IntersectionIterator nit = en.ibegin();

      double dtLocal = 0.0;
      double minvol = vol; 
      
      for (; nit != endnit; ++nit) {
	int twistSelf = twistUtil_.twistInSelf(nit); 
        FaceQuadratureType faceQuadInner(nit, quadOrd1_, twistSelf, 
                                         FaceQuadratureType::INSIDE);
	if (nit.neighbor()) {
	  EntityType& nb=*nit.outside();
          if (iset.index(nb) > iset.index(en)
              || nb.partitionType() == GhostEntity) {

            int twistNeighbor = twistUtil_.twistInNeighbor(nit);
            FaceQuadratureType faceQuadOuter(nit, quadOrd1_, twistNeighbor,
                                             FaceQuadratureType::OUTSIDE);
            
            caller_.setNeighbor(nb);
            LocalFunctionType updNeigh =dest_->localFunction(nb);

            const BaseFunctionSetType& bsetNeigh = 
              spc_.getBaseFunctionSet(nb);
  
	    double nbvol = volumeElement(nb, volQuad);
	    if (nbvol<minvol) minvol=nbvol;
            for (int l = 0; l < faceQuadInner.nop(); ++l) {
              double dtLocalS = 
                caller_.numericalFlux(nit, faceQuadInner, faceQuadOuter,
                                      l, valEn_, valNeigh_);
	      dtLocal += dtLocalS*faceQuadInner.weight(l);

              for (int i = 0; i < updEn.numDofs(); ++i) {
                updEn[i] -= 
                  bsetEn.evaluateSingle(i, faceQuadInner, l, valEn_)*
                  faceQuadInner.weight(l)/vol;
                updNeigh[i] += 
                  bsetNeigh.evaluateSingle(i, faceQuadOuter, l, valNeigh_)
                  *faceQuadOuter.weight(l)/nbvol;
              }
            }
                         
          } // end if ...
        } // end if neighbor

        if (nit.boundary()) {
          for (int l = 0; l < faceQuadInner.nop(); ++l) {
            double dtLocalS = 
              caller_.boundaryFlux(nit, faceQuadInner, l, source_);
	    dtLocal += dtLocalS*faceQuadInner.weight(l);
                    
            for (int i = 0; i < updEn.numDofs(); ++i) {
              updEn[i] -= bsetEn.evaluateSingle(i, faceQuadInner, l, source_)
                *faceQuadInner.weight(l)/vol;
            }
          }
        } // end if boundary
      }
      if (dtLocal>2.*std::numeric_limits<double>::min()) {
	dtLocal = minvol/dtLocal;
	if (dtLocal < dtMin_) dtMin_ = dtLocal;
      }
    }

  private:
    LocalDGPass();
    LocalDGPass(const LocalDGPass&);
    LocalDGPass& operator=(const LocalDGPass&);

  private:
    double volumeElement(const EntityType& en,
                         const VolumeQuadratureType& quad) const {
      double result = 0.0;
      for (int qp = 0; qp < quad.nop(); ++qp) {
        result += 
          quad.weight(qp) * en.geometry().integrationElement(quad.point(qp));
      }
      return result;
    }
    
  private:
    mutable DiscreteModelCallerType caller_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    DiscreteFunctionSpaceType& spc_;
    mutable double dtMin_;
  
    //! Some helper variables
    mutable JacobianRangeType fMat_;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType source_;
    mutable DomainType grads_;
    TimeProvider* time_;
    FieldVector<int, 0> diffVar_;

    TwistUtility<GridType> twistUtil_;

    int quadOrd0_,quadOrd1_;
  };
  
} // end namespace Dune

#endif
