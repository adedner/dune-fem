#ifndef DUNE_FEM_TIMEPROVIDER_HH
#define DUNE_FEM_TIMEPROVIDER_HH

//- system includes 
#include <limits>
#include <cassert>

//- Dune includes 
#include <dune/common/exceptions.hh>

#include <dune/fem/misc/commhelper.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/misc/femtuples.hh>

namespace Dune
{

  /** \class   TimeProviderBase
   *  \ingroup ODESolver
   *  \brief   general base for time providers
   *
   *  This class consists of the methods required for example in the
   *  ODE Solvers, e.g., provideTimeStepEstimate and 
   *  provideTimeStepUpperBound.
   *  InvalidateTimeStep can be used to mark this time step as invalid.
   *  Furthermore, method for accessing the simulation time, the
   *  time step counter and the time step size are provided.
   *
   *  The derived class TimeProvider provides the additional method
   *  required for implementing a time loop.
   *  
   */
  class TimeProviderBase : public AutoPersistentObject
  {
    typedef TimeProviderBase ThisType;

  protected:
    double time_;
    int timeStep_;
    double dt_;
    bool valid_;
    bool dtEstimateValid_;
    double dtEstimate_;
    double dtUpperBound_;

  public:
    inline TimeProviderBase ()
    : time_( Parameter :: getValue( "fem.timeprovider.starttime",
                                    (double)0.0 ) ),
      timeStep_( 0 ),
      dt_( 0.0 ),
      valid_( false ),
      dtEstimateValid_( false )
    {
      initTimeStepEstimate();
    }

    inline explicit TimeProviderBase ( const double startTime )
    : time_( startTime ),
      timeStep_( 0 ),
      dt_( 0.0 ),
      valid_( false ),
      dtEstimateValid_( false )
    {
      initTimeStepEstimate();
    }

    virtual ~TimeProviderBase() {}

    void backup() const 
    {
      Tuple<const double&,const int&,const double&,const bool&,const double&>
        values(time_,timeStep_,dt_,valid_,dtEstimate_);
      PersistenceManager::backupValue("timeprovider",values);
    }

    void restore() 
    {
      Tuple<double&,int&,double&,bool&,double&>
        values(time_,timeStep_,dt_,valid_,dtEstimate_);
      PersistenceManager::restoreValue("timeprovider",values);
      dtEstimateValid_ = true;
    }
    
  private:
    TimeProviderBase ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    /** \brief obtain the current time
     *
     *  \returns the current time
     */
    double time () const
    {
      return time_;
    }
    
    /** \brief obtain number of the current time step
     *
     *  \return the current time step counter
     */
    int timeStep () const
    {
      assert( timeStepValid() );
      return timeStep_;
    }
 
    /** \brief obtain the size of the current time step
     *
     *  \returns the size of the current time step
     */
    double deltaT () const
    {
      assert( timeStepValid() );
      return dt_;
    }

    /** \brief obtain current estimate on time step
     *
     *  \returns the current estimate for the time step
     */
    double timeStepEstimate () const
    {
      return dtEstimate_;
    }

    /** \brief set time step estimate to minimum of given value and
               internal time step estiamte 
         \param[in] dtEstimate time step size estimate 
    */
    void provideTimeStepEstimate ( const double dtEstimate )
    {
      dtEstimate_ = std::min( dtEstimate_, dtEstimate );
      dtEstimateValid_ = true;
    }
    /** \brief set upper bound for time step to minimum of given value and
               internal bound
         \param[in] upperBound time step size estimate 
    */
    void provideTimeStepUpperBound ( const double upperBound )
    {
      dtUpperBound_ = std::min( dtUpperBound_, upperBound );
      dtEstimateValid_ = true;
    }
    
    /** \brief count current time step a not valid */
    void invalidateTimeStep ()
    {
      valid_ = false;
    }

    /** \brief return if this time step should be used */
    bool timeStepValid () const
    {
      return valid_;
    }
   
  protected:
    void advance ()
    {
      if( timeStepValid() )
      {
        time_ += deltaT();
        ++timeStep_;
      }
    }

    void initTimeStepEstimate ()
    {
      dtEstimate_ = std::numeric_limits< double >::max();
      dtUpperBound_ = std::numeric_limits< double >::max();
      dtEstimateValid_ = false;
    }
  };



  /** 
   *  \ingroup ODESolver
   *  \brief   manager for global simulation time of time-dependent solutions
   *
   *  When calculating possibly multiple time-dependent solutions, it is often
   *  necessary to use the same time in all calculations. This means that we
   *  have to use the same time step for all our calculations. A TimeProvider
   *  keeps track of this information in a simple and unified way.
   *
   *  An example of a time loop could look as follows:
   *  \code
   *  // create time provider
   *  TimeProvider tp( startTime );
   *
   *  SpaceOperator spaceOperator;
   *  typedef SpaceOperator::DestinationType DestinationType;
   *  OdeSolver<DestinationType> odeSolver(spaceOperator,tp,order);
   *
   *  DestinationType U;
   *  initialize(U);
   *
   *  // set the initial time step estimate
   *  odeSolver.initialize( U );
   *
   *  // time loop
   *  for( tp.init(); tp.time() < endTime; tp.next() )
   *  {
   *    // do calculation
   *    odeSolver.solve(U);
   *  }
   *  \endcode
   *
   *  Within the time loop, both tp.time() and tp.deltaT() are fixed and cannot
   *  be altered and an the next time step should be fixed in the loop,
   *  e.g., in the method solve of the ode solver an upper estimate
   *  for the next time step is provided; if more than one time
   *  step restriction has to be imposed, the minimum is taken for
   *  the next time step.
   *  By calling the method provideTimeStepEstimate(maxDt) in the body of the
   *  loop an upper estimate for the next time step can be supplied;
   *  to fix the next time step (ignoring the estimates) an optinal
   *  argument can be passed to the next method on the
   *  Dune::TimeProvider.
   *
   *  Obviously, we need to provide an initial estimate. In the above example,
   *  this is done by the initialize method of the ODE solver. In tp.init(),
   *  the first time step (deltaT) is set based on the estimate and 
   *  this value can also be fixed independent of the estimate through
   *  an optional argument. The following loop would fix the time step
   *  to 1e-3
   *  \code
   *  for( tp.init(1e-3); tp.time() < endTime; tp.next(1e-3) )
   *  {
   *    // do calculation
   *    odeSolver.solve(U);
   *  }
   *  \endcode
   *
   *  In order to allow the user to incfluence the calculation of the next time
   *  step from the estimate, the time provider also maintains an additional
   *  factor (which is constant during the entire simulation). 
   *  Therefore the actual time step used, is calculated as follows:
   *  \f[
   *  \mathrm{deltaT} = \mathrm{factor} * \mathrm{timeStepEstimate}.
   *  \f]
   *  Therefore in the above example 1e-3 might not be the acctual
   *  time step depending on the value of the factor in the
   *  TimeProvider.
   *  The default value for this factor is equal to one but can be changed
   *  either during the construction of the Dune::TimeProvider or
   *  by using the parameter \c fem.timeprovider.factor.
   *  A further parameter read by the Dune::TimeProvider is
   *  fem.timeprovider.starttime defining the starting time of
   *    the simulation (default is zero).
   *
   *  The most general implementation is given in the class
   *  Dune::TimeProvider< CollectiveCommunication< C > >  which
   *  takes a Dune::CollectiveCommunication instance in the 
   *  constructor which is used in parallel computations is
   *  syncronize the time step. It defaults to 
   *  Dune::CollectiveCommHelperType :: defaultCommunication()
   *  and also works for seriell runs where the template argument
   *  does not have to be prescribed.
   *  If the communication manager from a given grid is to be used
   *  the class Dune::GridTimeProvider using the GridType as
   *  template argument can be used instead, with the same
   *  functionality.
   *
     \parametername \c fem.timeprovider.factor \n
                    multiplication factor to use for each time step;
                    defaults to 1.
     \parametername \c fem.timeprovider.starttime \n
                    time used for initializing the starting time
                    defaults to zero.
     \parametername \c fem.timeprovider.updatestep \n
                    only do the update of the time step size 
                    every 'updatestep' to avoid the 
                    expensive communication to achieve this 
                    (for testing only); 
                    defaults to 1 
   */
  template< class CommProvider = DefaultCollectiveCommunicationType >
  class TimeProvider {
  };
  
  /** \ingroup ODESolver
   *  \brief   the basic Dune::TimeProvider implementation.
   *
   *  This implementation of a timeprovider takes a CollectiveCommunicate 
   *  for parallel runs which default to a default communicator
   *  which also works for serial simulations.
   *
   */
  template< class C >
  class TimeProvider< CollectiveCommunication< C > >
  : public TimeProviderBase
  {
    typedef TimeProvider< CollectiveCommunication< C > > ThisType;
    typedef TimeProviderBase BaseType;

  public:
    typedef CollectiveCommunication< C > CollectiveCommunicationType;

  protected:
    typedef CollectiveCommunicationHelper< CollectiveCommunicationType >
      CollectiveCommHelperType;

    double getCflFactor() const 
    {
      return Parameter::getValidValue( "fem.timeprovider.factor", (double)1.0, ValidateGreater< double >( 0.0 ) );
    }
    
    int getUpdateStep () const 
    {
      return Parameter::getValidValue( "fem.timeprovider.updatestep", (int)1, ValidateGreater< int >( 0 ) );
    }
    
  public:
    /** \brief default constructor
     *
     *  \param[in]  comm  collective communication (optional)
     */
    explicit
    TimeProvider ( const CollectiveCommunicationType &comm
                     = CollectiveCommHelperType::defaultCommunication() )
    : BaseType(),
      comm_( comm ),
      cfl_( getCflFactor() ),
      updateStep_( getUpdateStep() ),
      counter_( updateStep_ )
    {}

    /** \brief constructor taking start time
     *
     *  \param[in]  startTime  initial time
     *  \param[in]  comm       collective communication (optional)
     
     */
    explicit
    TimeProvider ( const double startTime,
                   const CollectiveCommunicationType &comm
                     = CollectiveCommHelperType::defaultCommunication() )
    : BaseType( startTime ),
      comm_( comm ),
      cfl_( getCflFactor() ),
      updateStep_( getUpdateStep() ),
      counter_( updateStep_ )
    {}
    
    /** \brief constructor taking start time and CFL constant
     *
     *  \param[in]  startTime  initial time
     *  \param[in]  cfl        CFL constant
     *  \param[in]  comm       collective communication (optional)
     */
    TimeProvider ( const double startTime,
                   const double cfl,
                   const CollectiveCommunicationType &comm
                     = CollectiveCommHelperType :: defaultCommunication() )
    : BaseType( startTime ),
      comm_( comm ),
      cfl_( cfl ),
      updateStep_( 1 ),
      counter_( updateStep_ )
    {}

    virtual~TimeProvider() {}
    
  private:
    TimeProvider ( const ThisType & );
    ThisType &operator= ( const ThisType & );

  public:
    /** \brief init dt with time step estimate
     */
    void init ()
    {
      initTimeStep( dtEstimate_ );
    }

    /** \brief init dt with provided time step
     *
     *  \param[in]  timeStep  value of the first time step (is multiplied with
     *                        factor)
     */
    void init ( const double timeStep )
    {
      initTimeStep( timeStep );
    }
    
    /** \brief goto next time step
     *
     * Sets the size of the next time step to the current time step estimate
     * and sets the estimate to infinity.
     */
    void next () 
    {
      assert( this->dtEstimateValid_ );
      advance();
      initTimeStep( dtEstimate_ );
    }

    /** \brief goto next time step
     * 
     * Sets the size of the next time step to the provided time step value
     * and sets the estimate to infinity.
     * 
     *  \param[in]  timeStep  value of the next time step (is multiplied with
     *                        factor)
     */
    void next ( const double timeStep ) 
    {
      advance();
      initTimeStep(timeStep);
    }

    /** \brief  return the global factor number 
        \return time step factor 
    */
    double factor () const
    {
      return cfl_;
    }

  protected:
    using BaseType::advance;
    using BaseType::initTimeStepEstimate;

    void initTimeStep ( const double dtEstimate )
    {
      // increase counter 
      ++counter_ ;

      if( counter_ >= updateStep_ ) 
      {
        // set timestep estimate 
        dt_ = std::min(cfl_ * dtEstimate,dtUpperBound_);
        dt_ = comm_.min( dt_ );
        valid_ = (dt_ > 0.0);
        // reset counter 
        counter_ = 0;
      }

      initTimeStepEstimate();
    }

  public:
    /** \brief restore time and timestep from outside 
         (i.e. from former calculation)  
         \param[in] time new time 
         \param[in] timeStep new time step counter 
    */
    void restore ( const double time, const int timeStep )
    { 
      time_ = time; 
      timeStep_ = timeStep;
    }
    
    virtual void backup () const
    {
      BaseType::backup();
    }

    virtual void restore ()
    {
      BaseType::restore();
      const_cast< double & >( cfl_ ) = getCflFactor();
    }

  protected:
    using BaseType::dt_;
    using BaseType::dtEstimate_;
    using BaseType::dtUpperBound_;
    using BaseType::valid_;
    using BaseType::timeStep_;

    const CollectiveCommunicationType &comm_;
    const double cfl_;
    const int updateStep_;
    int counter_; 
  };

  /** \class   GridTimeProvider
   *  \ingroup ODESolver
   *  \brief   the same functionality as the Dune::TimeProvider.
   *
   *  This implementation of a timeprovider takes the CollectiveCommunicate 
   *  from a Dune::Grid instance.
   */
  template< class Grid >
  class GridTimeProvider
  : public TimeProvider< typename Grid::Traits::CollectiveCommunication >
  {
    typedef GridTimeProvider< Grid > ThisType;
    typedef TimeProvider< typename Grid::Traits::CollectiveCommunication > BaseType;

  public:
    typedef typename Grid::Traits::CollectiveCommunication CollectiveCommunicationType;

    explicit GridTimeProvider ( const Grid &grid )
    : BaseType( grid.comm() )
    {}

    GridTimeProvider ( const double startTime,
                       const Grid &grid )
    : BaseType( startTime, grid.comm() )
    {}
    
    GridTimeProvider ( const double startTime,
                       const double cfl,
                       const Grid &grid )
    : BaseType( startTime, cfl, grid.comm() )
    {}
    
    virtual ~GridTimeProvider() {}
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_TIMEPROVIDER_HH
