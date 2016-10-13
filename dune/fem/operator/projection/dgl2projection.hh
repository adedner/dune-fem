#ifndef DUNE_FEM_DGL2PROJECTION_HH
#define DUNE_FEM_DGL2PROJECTION_HH

#include <type_traits>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/operator/1order/localmassmatrix.hh>

namespace Dune
{

  namespace Fem
  {

    // DGL2ProjectionImpl
    // ------------------

    // implementation of L2 projection for discontinuous spaces
    class DGL2ProjectionImpl
    {
      template <int dummy, bool hasLocalFunction>
      struct ProjectChooser
      {
        template <class FunctionImp, class FunctionSpace>
        class FunctionAdapter
        {
          const FunctionImp& function_;
        public:
          typedef FunctionSpace FunctionSpaceType;
          typedef typename FunctionSpaceType :: RangeType RangeType;
          typedef typename FunctionSpaceType :: DomainType DomainType;
          FunctionAdapter(const FunctionImp& f) : function_(f) {}

          void evaluate(const DomainType& local,
                        RangeType& ret) const
          {
            function_.evaluate( local , ret );
          }
        };

        template <class FunctionImp, class DiscreteFunctionImp>
        static void project(const FunctionImp& f,
                            DiscreteFunctionImp& discFunc,
                            const int polOrd)
        {
          // some typedefs
          typedef typename DiscreteFunctionImp :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
          typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;
          typedef FunctionAdapter<FunctionImp, typename  DiscreteFunctionSpaceType :: FunctionSpaceType> FunctionAdapterType;
          // create function adapter in case of incorrect implementation
          FunctionAdapterType af( f );
          // create discrete function adapter
          GridFunctionAdapter< FunctionAdapterType, GridPartType> adapter(
              "DGL2projection::adapter" , f , discFunc.space().gridPart());
          DGL2ProjectionImpl::projectFunction(adapter, discFunc, polOrd);
        }
      };

      template <int dummy>
      struct ProjectChooser<dummy,true>
      {
        template <class FunctionImp, class DiscreteFunctionImp>
        static void project(const FunctionImp& f,
                            DiscreteFunctionImp& discFunc,
                            const int polOrd )
        {
          DGL2ProjectionImpl::projectFunction(f, discFunc, polOrd);
        }
      };

    public:
      /** /brief project function onto discrete discontinuous galerkin space
       *
       * \param f  function that is going to be projected
       * \param discFunc discrete function storing the result
       * \param quadOrd order of quadrature used (defaults to 2 * space.order())
       * \param communicate  restore integrity of data (defaults to true)
       */
      template <class FunctionImp, class DiscreteFunctionImp>
      static void project(const FunctionImp& f, DiscreteFunctionImp& discFunc,
                          const int quadOrd = -1, const bool communicate = true )
      {
        ProjectChooser<0, std::is_convertible<FunctionImp, HasLocalFunction>::value > :: project(f,discFunc,quadOrd);

        // do communication in parallel cases
        if( communicate )
          discFunc.communicate();
      }

      //! project function to discrete space
      template <class FunctionImp, class DiscreteFunctionImp>
      static void apply(const FunctionImp& f, DiscreteFunctionImp& discFunc )
      {
        project( f, discFunc );
      }

    protected:
      template <class FunctionImp, class DiscreteFunctionImp>
      static void projectFunction(const FunctionImp& func,
                                  DiscreteFunctionImp& discFunc,
                                  int polOrd = -1)
      {
        typedef typename DiscreteFunctionImp::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        typedef typename DiscreteFunctionImp::LocalFunctionType   LocalFuncType;
        typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
        typedef typename DiscreteFunctionSpaceType::RangeType     RangeType;

        typedef typename FunctionImp::LocalFunctionType LocalFType;

        const DiscreteFunctionSpaceType& space =  discFunc.space();

        // type of quadrature
        typedef CachingQuadrature<GridPartType,0> QuadratureType;
        // type of local mass matrix
        typedef LocalMassMatrix< DiscreteFunctionSpaceType, QuadratureType > LocalMassMatrixType;

        const int quadOrd = (polOrd == -1) ? (2 * space.order()) : polOrd;

        // create local mass matrix object
        LocalMassMatrixType massMatrix( space, quadOrd );

        // clear destination
        discFunc.clear();

        // extract type from grid part
        typedef typename GridPartType::template Codim<0>::GeometryType Geometry;

        // create storage for values
        std::vector< RangeType > values;

        for(const auto & en : space)
        {
          // get geometry
          const Geometry& geo = en.geometry();

          // get quadrature
          QuadratureType quad(en, quadOrd);

          // get local function of destination
          LocalFuncType lf = discFunc.localFunction(en);

          // get local function of argument
          const LocalFType f = func.localFunction(en);

          const int quadNop = quad.nop();

          // adjust size of values
          values.resize( quadNop );

          for(int qP = 0; qP < quadNop ; ++qP)
          {
            const double intel =
                 quad.weight(qP) * geo.integrationElement( quad.point(qP) );

            // evaluate function at quadrature point
            f.evaluate(quad[ qP ], values[ qP ] );

            // apply weight
            values[ qP ] *= intel;
          }

          // add values to local function
          lf.axpyQuadrature( quad, values );

          // apply inverse of mass matrix to local function
          massMatrix.applyInverse( en, lf );
        }
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DGL2PROJECTION_HH
