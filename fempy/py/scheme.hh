#ifndef DUNE_FEMPY_PY_SCHEME_HH
#define DUNE_FEMPY_PY_SCHEME_HH

#include <dune/fempy/pybind11/pybind11.hh>

#include <dune/common/typeutilities.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bcrsmatrix.hh>
#include <dune/python/istl/bcrsmatrix.hh>
#endif // #if HAVE_DUNE_ISTL

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/operator/matrix/colcompspmatrix.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/common/numpyvector.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/fempy/pybind11/pybind11.hh>

#if 0
     namespace pybind11
     {
       namespace detail
       {
          template <> class type_caster<_p_Mat>
          {
            public:
            PYBIND11_TYPE_CASTER(Mat, _("mat"));
            // Python to C++
            bool load(handle src, bool)
            {
              value = PyPetscMat_Get(src.ptr());
              return true;
            }
            static handle cast(Mat src, pybind11::return_value_policy policy, handle parent)
            {
               return pybind11::handle(PyPetscMat_New(src));
            }
            operator Mat() { return value; }
          };
        }
      }
#endif

namespace Dune
{

  namespace FemPy
  {

    // registerScheme
    // --------------

    namespace detail
    {

#if HAVE_DUNE_ISTL
      template< class B, class A >
      inline static const BCRSMatrix< B, A > &getBCRSMatrix ( const BCRSMatrix< B, A > &matrix ) noexcept
      {
        return matrix;
      }
#endif // #if HAVE_DUNE_ISTL



      // registerSchemeConstructor
      // -------------------------

      template< class Scheme, class... options >
      inline static auto registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_constructible< Scheme, const typename Scheme::DiscreteFunctionSpaceType &, const typename Scheme::ModelType & >::value >
      {
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Scheme::ModelType ModelType;

        using pybind11::operator""_a;

        cls.def( pybind11::init( [] ( Space &space, const ModelType &model ) {
            return new Scheme( space, model );
          } ), "space"_a, "model"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
        cls.def( pybind11::init( [] ( Space &space, const ModelType &model, const pybind11::dict &parameters ) {
            return new Scheme( space, model, pyParameter( parameters, std::make_shared< std::string >() ) );
          } ), "space"_a, "model"_a, "parameters"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeConstructor( cls, PriorityTag< 42 >() );
      }



      // registerSchemeAssemble
      // ----------------------
      // register assemble method if data method is available (and return value is registered)
#ifdef PETSC4PY_H // will be set it petsc4py.h was included (so import_petsc4py exists and the python module as well)
      template< class GF, class Scheme, class... options, std::enable_if_t<
            std::is_same< std::decay_t< decltype(std::declval< Scheme >().assemble( std::declval< const GF& >() )) >,
                  typename Scheme::LinearOperatorType>::value, int > _i = 0 >
      inline static auto registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 4 > )
        -> void_t< decltype( std::declval< const typename Scheme::LinearOperatorType & >().petscMatrix() ) >
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename DiscreteFunction::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;
        typedef typename Scheme::LinearOperatorType::MatrixType MatrixType;

        using pybind11::operator""_a;

        cls.def( "assemble", [] ( pybind11::handle obj, const GF &ubar ) {
            Scheme &self = obj.template cast<Scheme&>();

            if (import_petsc4py() != 0)
            {                           \
              std::cout << "ERROR: could not import petsc4py\n";
              throw std::runtime_error("Error during import of petsc4py");
            }

            Mat mat = self.assemble( ubar ).petscMatrix();
            pybind11::handle petsc_mat(PyPetscMat_New(mat));
            return petsc_mat;
          }, "ubar"_a );
      }
#endif
#if HAVE_DUNE_ISTL
      template< class GF, class Scheme, class... options, std::enable_if_t<
            std::is_same< std::decay_t< decltype(std::declval< Scheme >().assemble( std::declval< const GF& >() )) >,
                  typename Scheme::LinearOperatorType>::value, int > i = 0 >
      inline static auto registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 3 > )
        -> void_t< decltype( getBCRSMatrix( std::declval< const typename Scheme::LinearOperatorType & >().matrix() ) ) >
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename DiscreteFunction::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;

        typedef std::decay_t< decltype( getBCRSMatrix( std::declval< const typename Scheme::LinearOperatorType & >().matrix() ) ) > BCRSMatrix;
        if( !pybind11::already_registered< BCRSMatrix >() )
          Python::registerBCRSMatrix< BCRSMatrix >( cls );

        using pybind11::operator""_a;

        cls.def( "assemble", [] ( Scheme &self, const GF &ubar ) {
            return getBCRSMatrix( self.assemble( ubar ).matrix() );
          }, pybind11::return_value_policy::reference_internal, "ubar"_a );
      }
#endif // #if HAVE_DUNE_ISTL

      template< class GF, class Scheme, class... options, std::enable_if_t<
            std::is_same< std::decay_t< decltype(std::declval< Scheme >().assemble( std::declval< const GF& >() )) >,
                  typename Scheme::LinearOperatorType>::value, int > i = 0 >
      inline static auto registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 2 > )
        -> void_t< decltype( std::declval< const typename Scheme::LinearOperatorType & >().matrix().data() ) >
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename DiscreteFunction::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;

        using pybind11::operator""_a;

        cls.def( "assemble", [] ( Scheme &self, const GF &ubar ) {
            return self.assemble( ubar ).matrix().data();
          }, pybind11::return_value_policy::reference_internal, "ubar"_a );
      }

      template< class GF, class Scheme, class... options, std::enable_if_t<
            std::is_same< std::decay_t< decltype(std::declval< Scheme >().assemble( std::declval< const GF& >() )) >,
                  typename Scheme::LinearOperatorType>::value, int > _i = 0 >
      inline static auto registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< const typename Scheme::LinearOperatorType & >().matrix() ) >
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename DiscreteFunction::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;
        typedef typename Scheme::LinearOperatorType::MatrixType MatrixType;

        using pybind11::operator""_a;

        cls.def( "assemble", [] ( pybind11::handle obj, const GF &ubar ) {
            Scheme &self = obj.template cast<Scheme&>();

            auto& mat = self.assemble( ubar ).matrix();

            pybind11::array_t<size_t> outerIndices(mat.rows() + 1);

            size_t nnz = 0;
            for(size_t i=0;i<mat.rows();++i)
              nnz += mat.numNonZeros(i);

            pybind11::array_t<size_t> innerIndices(nnz);
            pybind11::array_t<double> data(nnz);

            size_t fill = 0;
            outerIndices.mutable_at(0) = 0;
            for(size_t i=0;i<mat.rows();++i)
            {
              size_t count = i*mat.numNonZeros();
              for(size_t j=0;j<mat.numNonZeros();++j,++count)
              {
                const auto pairIdx = mat.realValue(count);
                if (pairIdx.second < mat.cols())
                {
                  innerIndices.mutable_at(fill) = pairIdx.second;
                  data.mutable_at(fill) = pairIdx.first;
                  ++fill;
                }
                else break;
              }
              outerIndices.mutable_at(i+1) = fill;
            }
            pybind11::object matrix_type = pybind11::module::import("scipy.sparse").attr("csr_matrix");
            pybind11::object scipy_mat = matrix_type(
                std::make_tuple(data, innerIndices, outerIndices),
                std::make_pair(mat.rows(), mat.cols())
            );
            return scipy_mat;
          }, "ubar"_a );
      }
      template< class GF, class Scheme, class... options >
      inline static void registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {
      }

      template< class Scheme, class... options >
      inline static void registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls )
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename DiscreteFunction::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;

        registerSchemeAssemble< DiscreteFunction >( cls, PriorityTag< 42 >() );
        registerSchemeAssemble< VirtualizedGridFunction< GridPart, RangeType > >( cls, PriorityTag< 42 >() );
      }



      // registerSchemeGeneralCall
      // -------------------------

      template< class Scheme, class... options >
      inline static auto registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< Scheme & >()(
                     std::declval< const VirtualizedGridFunction< typename Scheme::GridPartType, typename Scheme::DiscreteFunctionSpaceType::RangeType > & >(),
                     std::declval< typename Scheme::DiscreteFunctionType & >()
                   ) ) >
      {
        typedef typename Scheme::DiscreteFunctionSpaceType::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        cls.def( "__call__", [] ( Scheme &self, const VirtualizedGridFunction< GridPart, RangeType > &arg, DiscreteFunction &dest ) { self( arg, dest ); } );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeGeneralCall( cls, PriorityTag< 42 >() );
      }



      // registerSchemeModel
      // -------------------

      template< class Scheme, class... options >
      inline static auto registerSchemeModel ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< Scheme >().model() ) >
      {
        cls.def_property_readonly( "model", &Scheme::model );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeModel ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeModel ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeModel( cls, PriorityTag< 42 >() );
      }



      // registerScheme
      // --------------

      template< class Scheme, class... options >
      inline static void registerScheme ( pybind11::module module, pybind11::class_< Scheme, options... > cls )
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;

        using pybind11::operator""_a;

        registerSchemeConstructor( cls );

        cls.def( "_solve", [] ( Scheme &self, const DiscreteFunction &rhs, DiscreteFunction &solution ) {
            auto info = self.solve( rhs, solution );
            // needs pybind 1.9: return pybind11::dict("converged"_a=info.converged, "iterations"_a=info.nonlinearIterations, "linear_iterations"_a=info.linearIterations);
            return std::map<std::string,std::string> {
                {"converged",std::to_string(info.converged)},
                {"iterations",std::to_string(info.nonlinearIterations)},
                {"linear_iterations",std::to_string(info.linearIterations)}
              };
          } );
        cls.def( "_solve", [] ( Scheme &self, DiscreteFunction &solution ) {
            auto info = self.solve( solution );
            // needs pybind 1.9: return pybind11::dict("converged"_a=info.converged, "iterations"_a=info.nonlinearIterations, "linear_iterations"_a=info.linearIterations);
            return std::map<std::string,std::string> {
                {"converged",std::to_string(info.converged)},
                {"iterations",std::to_string(info.nonlinearIterations)},
                {"linear_iterations",std::to_string(info.linearIterations)}
              };
          } );
        cls.def( "__call__", [] ( Scheme &self, const DiscreteFunction &arg, DiscreteFunction &dest) { self( arg, dest ); } );
        registerSchemeGeneralCall( cls );

        cls.def_property_readonly( "dimRange", [] ( Scheme & ) -> int { return DiscreteFunction::FunctionSpaceType::dimRange; } );
        cls.def_property_readonly( "space", [] ( pybind11::object self ) { return detail::getSpace( self.cast< const Scheme & >(), self ); } );
        registerSchemeModel( cls );

        registerSchemeAssemble( cls );

        cls.def( "constraint", [] ( Scheme &self, DiscreteFunction &u) { self.constraint( u ); } );

        cls.def( "mark", [] ( Scheme &self, const DiscreteFunction &solution, double tolerance ) {
            double est = self.estimate( solution );
            return std::make_tuple( est, self.mark( tolerance ) );
          } );
      }

    } // namespace detail

    template< class Scheme, class... options >
    inline static void registerScheme ( pybind11::module module, pybind11::class_< Scheme, options... > cls )
    {
      detail::registerScheme( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SCHEME_HH
