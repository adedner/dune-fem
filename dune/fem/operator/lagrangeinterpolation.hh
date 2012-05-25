#ifndef DUNE_FEM_LAGRANGEINTERPOLATION_HH
#define DUNE_FEM_LAGRANGEINTERPOLATION_HH

#include <dune/common/typetraits.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>

namespace Dune 
{

  /** \class LagrangeInterpolation
   *  \brief Generates the Lagrange Interpolation of an analytic function
   */
  template< class DiscreteFunction >
  class LagrangeInterpolation
  {
    typedef LagrangeInterpolation< DiscreteFunction > ThisType;

  public:
    //! type of discrete functions
    typedef DiscreteFunction DiscreteFunctionType;

    //! type of discrete function space
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;
    //! type of local functions
    typedef typename DiscreteFunctionType::LocalFunctionType
      LocalFunctionType;

    //! type of grid partition
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    //! type of Lagrange point set
    typedef typename DiscreteFunctionSpaceType::LagrangePointSetType
      LagrangePointSetType;
    //! type of vectors in function's domain
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    //! type of vectors in function's range
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

  public:
    /** interpolate an analytical function into a Lagrange discrete function
     *
     *  This Method evaluates the given function (which can be evaluated
     *  globally) at the Lagrange points and writes the values into a discrete
     *  function.
     *
     *  \param[in] function function to interpolate
     *
     *  \param[out] discreteFunction discrete function to receive the
     *              interpolation
     */
    template< class Function >
    static void interpolateFunction ( const Function &function,
                                      DiscreteFunctionType &discreteFunction );

    //! \copydoc interpolateFunction 
    //- make interface equal to DGL2Projection 
    template< class Function >
    static void apply ( const Function &function,
                        DiscreteFunctionType &discreteFunction )
    {
      interpolateFunction( function, discreteFunction );
    }

  private:
    template< class Function, bool hasLocalFunction >
    struct CallInterpolateDiscreteFunction;

  protected:
    /** interpolate a discrete function into a Lagrange discrete function
     *
     *  This Method evaluates the given function (which can be evaluated
     *  locally at the Lagrange points and writes the values into a discrete
     *  function.
     *
     *  \param[in] function function to interpolate
     *
     *  \param[out] discreteFunction discrete function to receive the
     *              interpolation
     */
    template< class Function >
    static void
    interpolateDiscreteFunction ( const Function &function,
                                  DiscreteFunctionType &discreteFunction );
  };



  template< class DiscreteFunction >
  template< class Function >
  inline void LagrangeInterpolation< DiscreteFunction >
    ::interpolateFunction ( const Function &function,
                             DiscreteFunctionType &discreteFunction )
  {
    const bool hasLocalFunction = Conversion< Function, HasLocalFunction >::exists;
    CallInterpolateDiscreteFunction< Function, hasLocalFunction >::call( function, discreteFunction );
  }



  template< class DiscreteFunction >
  template< class Function >
  struct LagrangeInterpolation< DiscreteFunction >
    ::CallInterpolateDiscreteFunction< Function, true >
  {
    static void call( const Function &function,
                      DiscreteFunction &discreteFunction )
    {
      interpolateDiscreteFunction( function, discreteFunction );
    }
  };

  template< class DiscreteFunction >
  template< class Function >
  struct LagrangeInterpolation< DiscreteFunction >
    ::CallInterpolateDiscreteFunction< Function, false >
  {
    static void call ( const Function &function,
                       DiscreteFunction &discreteFunction )
    {
      typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
      typedef GridFunctionAdapter< Function, GridPartType > GridFunctionAdapterType;

      const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();
      GridFunctionAdapterType dfAdapter( "function", function, dfSpace.gridPart() );
      interpolateDiscreteFunction( dfAdapter, discreteFunction );
    }
  };


  
  template< class DiscreteFunction >
  template< class Function >
  inline void LagrangeInterpolation< DiscreteFunction >
    ::interpolateDiscreteFunction ( const Function &function,
                                    DiscreteFunctionType &discreteFunction )
  {
    typedef typename DiscreteFunctionType::DofType DofType;
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
    typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
    static const int dimRange = DiscreteFunctionSpaceType::dimRange;

    typedef typename Function::LocalFunctionType FunctionLocalFunctionType;

    // set all DoFs to infinity
    const DofIteratorType dend = discreteFunction.dend();
    for( DofIteratorType dit = discreteFunction.dbegin(); dit != dend; ++dit )
      *dit = std::numeric_limits< DofType >::infinity();

    const DiscreteFunctionSpaceType &dfSpace = discreteFunction.space();

    IteratorType endit = dfSpace.end();
    for( IteratorType it = dfSpace.begin(); it != endit; ++it )
    {
      const LagrangePointSetType &lagrangePointSet
        = dfSpace.lagrangePointSet( *it );

      FunctionLocalFunctionType f_local = function.localFunction( *it );
      LocalFunctionType df_local = discreteFunction.localFunction( *it );

      // assume point based local dofs 
      const int nop = lagrangePointSet.nop();
      int k = 0;
      for( int qp = 0; qp < nop; ++qp )
      {
        // if the first DoF for this point is already valid, continue
        if( df_local[ k ] == std::numeric_limits< DofType >::infinity() )
        {
          // evaluate the function in the Lagrange point
          RangeType phi;
          f_local.evaluate( lagrangePointSet[ qp ], phi );

          // assign the appropriate values to the DoFs
          for( int i = 0; i < dimRange; ++i, ++k )
            df_local[ k ] = phi[ i ];
        }
        else
          k += dimRange;
      }
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_FEM_LAGRANGEINTERPOLATION_HH
