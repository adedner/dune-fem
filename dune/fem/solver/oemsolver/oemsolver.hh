#ifndef DUNE_FEM_OEMSOLVER_HH
#define DUNE_FEM_OEMSOLVER_HH

//- system includes
#include <limits>
#include <type_traits>
#include <utility>

//- Dune includes
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>

//- local includes
#include "preconditioning.hh"

#include <dune/fem/solver/pardg.hh>

// include BLAS  implementation
#include "cblas.h"


namespace Dune
{

  namespace Fem
  {

    /** \class   OEMMatrix
     *  \ingroup OEMSolver
     *  \brief   interface for matrices to be used with OEM sovlers
     */
    struct OEMMatrix
    {
      virtual ~OEMMatrix () {}

      /** \brief evaluate matrix vector multiplication
       *
       *  \param[in]   u  vector to multiply the matrix with
       *  \param[out]  w  vector to store the result in
       */
      virtual void multOEM ( const double *u, double *w ) const = 0;

      /** \brief evaluate scalar product
       *
       *  \param[in]   u  first argument of scalar product
       *  \param[in]   v  second argument of scalar product
       */
      virtual double ddotOEM ( const double *u, const double *v ) const = 0;

      // SparseRowMatrixObject does not implement this method
#if 0
      /** \brief apply preconditioner
       *
       *  \param[in]   u  vector to apply preconditioner to
       *  \param[out]  w  vector to store the result in
       */
      virtual void precondition ( const double *u, double *w ) const = 0;
#endif
    };

    } // end namespace Fem

  } // end namespace Dune



  namespace OEMSolver
  {

  //////////////////////////////////////////////////////////
  //
  // Operator Interface to use linear solvers from pardg
  //
  //////////////////////////////////////////////////////////
  template <class OperatorImp>
  class SolverInterfaceImpl
#ifdef USE_PARDG_ODE_SOLVER
  : public PARDG::Function
#endif
  {
    const OperatorImp & op_;
    int size_;
  public:
    SolverInterfaceImpl(const OperatorImp & op, int size = 0)
      : op_(op), size_(size)
    {}

    void setSize( int size ) { size_ = size; }

    void operator () (const double *arg, double * dest, int i = 0 )
    {
      op_.multOEM(arg,dest);
    }

    void mult(const double *arg, double * dest) const
    {
      op_.multOEM(arg,dest);
    }

    int dim_of_argument(int i = 0) const
    {
      assert( i == 0 );
      return size_;
    }
    int dim_of_value(int i = 0) const
    {
      assert( i == 0 );
      return size_;
    }
  };

  //////////////////////////////////////////////////////////
  //
  // Preconditioner Interface to use linear solvers from pardg
  //
  //////////////////////////////////////////////////////////
  template <class PreconditionerImp>
  class PreconditionerImpl
#ifdef USE_PARDG_ODE_SOLVER
  : public PARDG::Function
#endif
  {
    const PreconditionerImp& pre_;
    int size_;
  public:
    PreconditionerImpl(const PreconditionerImp& pre, int size = 0)
      : pre_(pre), size_(size)
    {}

    void setSize( int size ) { size_ = size; }

    void operator () (const double *arg, double * dest, int i = 0 )
    {
      pre_.precondition(arg,dest);
    }

    void mult(const double *arg, double * dest) const
    {
      pre_.precondition(arg,dest);
    }

    int dim_of_argument(int i = 0) const
    {
      assert( i == 0 );
      return size_;
    }
    int dim_of_value(int i = 0) const
    {
      assert( i == 0 );
      return size_;
    }
  };

  // use cblas implementations
  using namespace DuneCBlas;

  using DuneCBlas :: daxpy;
  using DuneCBlas :: dcopy;
  using DuneCBlas :: ddot;
  using DuneCBlas :: dnrm2;
  using DuneCBlas :: dscal;

  //! this method is called from all solvers and is only a wrapper
  //! this method is mainly from SparseRowMatrix
  template <class MatrixImp, class VectorType>
  void mult(const MatrixImp & m, const VectorType * x, VectorType * ret)
  {
    // call multOEM of the matrix
    m.multOEM(x,ret);
  }

  //! mult method when given pre conditioning matrix
  template <class Matrix , class PC_Matrix , bool >
  struct Mult
  {
    static inline double ddot( const Matrix& A,
                               const double *x,
                               const double *y)
    {
      return A.ddotOEM(x,y);
    }

    typedef bool mult_t(const Matrix &A,
                        const PC_Matrix & C,
                        const double *arg,
                        double *dest ,
                        double * tmp);

    static void back_solve(const int size,
          const PC_Matrix & C, double* solution, double* tmp)
    {
      assert( tmp );
      if( C.rightPrecondition() )
      {
        C.precondition(solution,tmp);
        // copy modified solution
        std::memcpy(solution,tmp, size * sizeof(double));
      }
    }

    static bool mult_pc (const Matrix &A, const PC_Matrix & C,
          const double *arg, double *dest , double * tmp)
    {
      assert( tmp );

      bool rightPreCon = C.rightPrecondition();
      // check type of preconditioning
      if( rightPreCon )
      {
        // call precondition of Matrix PC
        C.precondition(arg,tmp);

        // call mult of Matrix A
        mult(A,tmp,dest);
      }
      else
      {
        // call mult of Matrix A
        mult(A,arg,tmp);

        // call precondition of Matrix PC
        C.precondition(tmp,dest);
      }
      return rightPreCon ;
    }
  };

  //! mult method when no pre conditioning matrix
  template <class Matrix>
  struct Mult<Matrix,Matrix,false>
  {
    static inline double ddot( const Matrix& A,
                               const double *x,
                               const double *y)
    {
      return A.ddotOEM(x,y);
    }

    typedef bool mult_t(const Matrix &A,
                        const Matrix &C,
                        const double *arg,
                        double *dest ,
                        double * tmp);

    static void back_solve(const int size,
          const Matrix & C, double* solution, double* tmp)
    {
      // do nothing here
    }

    static bool mult_pc(const Matrix &A, const Matrix & C, const double *arg ,
                        double *dest , double * tmp)
    {
      // tmp has to be 0
      assert( tmp == 0 );
      // C is just a fake
      assert( &A == &C );

      // call mult of Matrix A
      mult(A,arg,dest);

      return true;
    }
  };

#define USE_MEMPROVIDER
#include "bicgstab.h"
#include "cghs.h"
#include "gmres.h"
#include "bicgsq.h"
#undef USE_MEMPROVIDER


  //! fake conditioner which just is id for internal parts of vector and zero
  //! for other parts, needed by parallel gmres
  class FakeConditioner
  {
    // size of vectors
    const int size_;

    // indices of external values
    std::vector<int> indices_;
  public:
    // use with care, not sure that working correctly already
    template <class SolverOperatorImp>
    FakeConditioner(int size, SolverOperatorImp& op) : size_(size)
    {
      assert( size_ > 0 );

      double * diag  = new double [size_];
      double * tmp   = new double [size_];

      assert( diag );
      assert( tmp );
      for(int i=0; i<size_; ++i) tmp[i] = i;
      op(tmp,diag);

      int newSize = (int) 0.25 * size_;
      indices_.reserve( newSize );
      indices_.resize( 0 );
      // now diag contains only non-zeros for all internal entries
      // these are set to 1.0 to be the id mapping
      for(int i=0; i<size_; ++i)
      {
        if( ! (std::abs (diag[i]) > 0.0) )
        {
          indices_.push_back( i );
        }
      }

      delete [] diag;
      delete [] tmp;
    }

    bool rightPrecondition() const { return false; }

    //! only keep internal parts of arg
    void precondition(const double * arg, double * dest) const
    {
      multOEM(arg,dest);
    }

    //! only keep internal parts of arg
    void multOEM(const double * arg, double * dest) const
    {
      std::memcpy( dest, arg , size_ * sizeof(double) );

      const int s = indices_.size();
      for(int i=0; i<s; ++i)
      {
        dest[indices_[i]] = 0.0;
      }
    }
  };

  } // end namespace OEMSolver

  namespace Dune
  {
    namespace Fem
    {
    /** @addtogroup OEMSolver

        In this section implementations of Orthogonal Error Methods (OEM) for solving linear
        systems of the from \f$A x = b\f$, where \f$A\f$ is a Mapping or
        Operator and \f$x\f$ and \f$b\f$ are discrete functions
        (see DiscreteFunctionInterface) can be found.

        @{
     **/

    /** \brief OEM-CG scheme after Hestenes and Stiefel */
    template <class DiscreteFunctionType, class OpType >
    class OEMCGOp : public Operator<
                DiscreteFunctionType,DiscreteFunctionType>
    {
    public:
      typedef OpType OperatorType;

    private:
      // no const reference, we make const later
      OperatorType *op_ = nullptr;
      typename DiscreteFunctionType::RangeFieldType epsilon_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;

      typedef std::pair < int , double > ReturnValueType;

      template <class OperatorImp, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class DiscreteFunctionImp>
        static ReturnValueType call(OperatorImp & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest,
                         double eps, int maxIter, bool verbose)
        {
          // use communication class of grid
          // see dune-common/common/collectivecommunication.hh
          // for interface
          int size = arg.space().size();

          if(op.hasPreconditionMatrix())
          {
            return OEMSolver::cghs(arg.space().gridPart().comm(),
                       size,op.systemMatrix(),op.preconditionMatrix(),
                       arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
          else
          {
            return OEMSolver::cghs(arg.space().gridPart().comm(),
                      size,op.systemMatrix(),
                      arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
        }
      };

      //! without any preconditioning
      template <class OperatorImp>
      struct SolverCaller<OperatorImp,false>
      {
        template <class DiscreteFunctionImp>
        static ReturnValueType call(OperatorImp & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest,
                         double eps, int maxIter, bool verbose)
        {
          // use communication class of grid
          // see dune-common/common/collectivecommunication.hh
          // for interface
          int size = arg.space().size();
          return OEMSolver::cghs(arg.space().gridPart().comm(),
                    size,op.systemMatrix(),
                    arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
        }
      };

    public:

      /** \brief constructor of OEM-CG
          \param[in] redEps realative tolerance for residual
          \param[in] absLimit absolut solving tolerance for residual
          \param[in] maxIter maximal number of iterations performed
          \param[in] verbose verbosity
      */
      OEMCGOp ( double redEps, double absLimit, int maxIter, bool verbose,
                const ParameterReader &parameter = Parameter::container() )
      : epsilon_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 )
      {}

      OEMCGOp ( double redEps, double absLimit, int maxIter,
                const ParameterReader &parameter = Parameter::container() )
      : OEMCGOp( redEps, absLimit, maxIter, parameter.getValue< bool >( "fem.solver.verbose", false ), parameter )
      {}

      OEMCGOp ( double redEps, double absLimit,
                const ParameterReader &parameter = Parameter::container() )
      : OEMCGOp( redEps, absLimit, std::numeric_limits< int >::max(),
                 parameter.getValue< bool >( "fem.solver.verbose", false ), parameter )
      {}

      /** \brief constructor of OEM-CG
          \param[in] op Operator to invert
          \param[in] redEps realative tolerance for residual
          \param[in] absLimit absolut solving tolerance for residual
          \param[in] maxIter maximal number of iterations performed
          \param[in] verbose verbosity
      */
      OEMCGOp ( OperatorType &op,
                double redEps, double absLimit, int maxIter, bool verbose,
                const ParameterReader &parameter = Parameter::container() )
      : OEMCGOp( redEps, absLimit, maxIter_, verbose, parameter )
      {
        bind( op );
      }

      OEMCGOp ( OperatorType &op,
                double redEps, double absLimit, int maxIter,
                const ParameterReader &parameter = Parameter::container() )
      : OEMCGOp( redEps, absLimit, maxIter_, parameter )
      {
        bind( op );
      }

      OEMCGOp ( OperatorType &op,
                double redEps, double absLimit,
                const ParameterReader &parameter = Parameter::container() )
      : OEMCGOp( redEps, absLimit, parameter )
      {
        bind( op );
      }

      void bind ( OperatorType& op )
      {
        op_ = &op;
      }

      void unbind ()
      {
        op_ = nullptr;
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {
      }

      int iterations ()  const
      {
        return iterations_;
      }

      void setMaxIterations ( int maxIter )
      {
        maxIter_ = maxIter;
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        assert( op_ );
        // prepare operator
        prepare ( arg, dest );

        ReturnValueType val =
          SolverCaller<OperatorType,
                       // check wheter operator has precondition methods
                       // to enable preconditioning derive your operator from
                       // OEMSolver::PreconditionInterface
                       std::is_convertible<OperatorType, OEMSolver::PreconditionInterface > ::value >::
                         // call solver, see above
                         call(*op_,arg,dest,epsilon_,maxIter_,verbose_);

        iterations_ = val.first;

        if(arg.space().gridPart().comm().rank() == 0)
        {
          std::cout << "OEM-CG: " << val.first << " iterations! Error: " << val.second << "\n";
        }

        // finalize operator
        finalize ();
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }
    };

    /** \brief BiCG-stab solver */
    template <class DiscreteFunctionType, class OpType>
    class OEMBICGSTABOp : public Operator<
                DiscreteFunctionType,DiscreteFunctionType>
    {
    public:
      typedef OpType OperatorType;

    private:
      // no const reference, we make const later
      OperatorType *op_ = nullptr;
      typename DiscreteFunctionType::RangeFieldType epsilon_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;

      typedef std::pair < int , double > ReturnValueType;

      template <class OperatorImp, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class DiscreteFunctionImp>
        static ReturnValueType call(OperatorImp & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest,
                         double eps, int maxIter, bool verbose)
        {
          int size = arg.space().size();
          if(op.hasPreconditionMatrix())
          {
            return OEMSolver::bicgstab(arg.space().gridPart().comm(),
                      size,op.systemMatrix(),op.preconditionMatrix(),
                      arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
          else
          {
            return OEMSolver::bicgstab(arg.space().gridPart().comm(),
                      size,op.systemMatrix(),
                      arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
        }
      };

      //! without any preconditioning
      template <class OperatorImp>
      struct SolverCaller<OperatorImp,false>
      {
        template <class DiscreteFunctionImp>
        static ReturnValueType call(OperatorImp & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest,
                         double eps, int maxIter, bool verbose)
        {
          int size = arg.space().size();
          return OEMSolver::bicgstab(arg.space().gridPart().comm(),
                    size,op.systemMatrix(),
                    arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
        }
      };

    public:
      /** \brief constructor of OEM-BiCG-stab
          \param[in] redEps realative tolerance for residual
          \param[in] absLimit absolut solving tolerance for residual
          \param[in] maxIter maximal number of iterations performed
          \param[in] verbose verbosity
      */
      OEMBICGSTABOp ( double redEps, double absLimit, int maxIter, bool verbose,
                      const ParameterReader &parameter = Parameter::container() )
      : epsilon_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 )
      {
      }

      OEMBICGSTABOp ( double redEps, double absLimit, int maxIter,
                      const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSTABOp( redEps, absLimit, maxIter,
                       parameter.getValue< bool >( "fem.solver.verbose", false ), parameter )
      {
      }

      OEMBICGSTABOp ( double redEps, double absLimit,
                      const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSTABOp( redEps, absLimit, std::numeric_limits< int >::max(),
                       parameter.getValue< bool >( "fem.solver.verbose", false ), parameter )
      {
      }

      /** \brief constructor of OEM-BiCG-stab
          \param[in] op Operator to invert
          \param[in] redEps realative tolerance for residual
          \param[in] absLimit absolut solving tolerance for residual
          \param[in] maxIter maximal number of iterations performed
          \param[in] verbose verbosity
      */
      OEMBICGSTABOp ( OperatorType& op,
                      double redEps, double absLimit, int maxIter, bool verbose,
                      const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSTABOp( redEps, absLimit, maxIter, verbose, parameter )
      {
        bind( op );
      }

      OEMBICGSTABOp ( OperatorType &op,
                      double redEps, double absLimit, int maxIter,
                      const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSTABOp( redEps, absLimit, maxIter, parameter )
      {
        bind( op );
      }

      OEMBICGSTABOp ( OperatorType &op,
                      double redEps, double absLimit,
                      const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSTABOp( redEps, absLimit, parameter )
      {
        bind( op );
      }

      void bind ( OperatorType &op )
      {
        op_ = &op;
      }

      void unbind ()
      {
        op_ = nullptr;
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {
      }

      int iterations () const
      {
        return iterations_;
      }

      void setMaxIterations ( int maxIter )
      {
        maxIter_ = maxIter;
      }


      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        assert( op_ );
        // prepare operator
        prepare ( arg, dest );

        ReturnValueType val =
          SolverCaller<OperatorType,
                       // check wheter operator has precondition methods
                       // to enable preconditioning derive your operator from
                       // OEMSolver::PreconditionInterface
                       std::is_convertible<OperatorType, OEMSolver::PreconditionInterface > ::value >::
                         // call solver, see above
                         call(*op_,arg,dest,epsilon_,maxIter_,verbose_);

        iterations_ = val.first;

        if(arg.space().gridPart().comm().rank() == 0)
        {
          std::cout << "OEM-BICGstab: " << val.first << " iterations! Error: " << val.second << "\n";
        }

        // finalize operator
        finalize ();
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }

    };

    ////////////////////////////////
    // BICG SQ scheme
    ////////////////////////////////
    /** \brief BiCG-SQ method */
    template <class DiscreteFunctionType, class OpType>
    class OEMBICGSQOp : public Operator<
                DiscreteFunctionType,DiscreteFunctionType>
    {
    public:
      typedef OpType OperatorType;

    private:
      // no const reference, we make const later
      OperatorType *op_ = nullptr;
      typename DiscreteFunctionType::RangeFieldType epsilon_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;

    public:
      /** \brief constructor of OEM-BiCG-SQ
          \param[in] op Operator to invert
          \param[in] redEps realative tolerance for residual
          \param[in] absLimit absolut solving tolerance for residual
          \param[in] maxIter maximal number of iterations performed
          \param[in] verbose verbosity
      */
      OEMBICGSQOp ( double redEps, double absLimit, int maxIter, bool verbose,
                    const ParameterReader &parameter = Parameter::container() )
      : epsilon_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 )
      {
      }

      OEMBICGSQOp ( double redEps, double absLimit,
                    const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSQOp( redEps, absLimit, std::numeric_limits< int >::max(),
                     parameter.getValue< bool >( "fem.solver.verbose", false ), parameter )
      {
      }

      OEMBICGSQOp ( double redEps, double absLimit, int maxIter,
                    const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSQOp( redEps, absLimit, maxIter,
                     parameter.getValue< bool >( "fem.solver.verbose", false ), parameter )
      {
      }

      /** \brief constructor of OEM-BiCG-SQ
          \param[in] op Operator to invert
          \param[in] redEps realative tolerance for residual
          \param[in] absLimit absolut solving tolerance for residual
          \param[in] maxIter maximal number of iterations performed
          \param[in] verbose verbosity
      */
      OEMBICGSQOp ( OperatorType &op,
                    double redEps, double absLimit, int maxIter, bool verbose,
                    const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSQOp( redEps, absLimit, maxIter, verbose, parameter )
      {
        bind( op );
      }

      OEMBICGSQOp ( OperatorType &op,
                    double redEps, double absLimit,
                    const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSQOp( redEps, absLimit, parameter )
      {
        bind( op );
      }

      OEMBICGSQOp ( OperatorType &op,
                    double redEps, double absLimit, int maxIter,
                    const ParameterReader &parameter = Parameter::container() )
      : OEMBICGSQOp( redEps, absLimit, maxIter, parameter )
      {
        bind( op );
      }

      void bind ( OperatorType &op )
      {
        op_ = &op;
      }

      void unbind ()
      {
        op_ = nullptr;
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {
      }

      int iterations () const
      {
        return iterations_;
      }

      void setMaxIterations ( int maxIter )
      {
        maxIter_ = maxIter;
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        assert( op_ );
        // prepare operator
        prepare ( arg, dest );

        int size = arg.space().size();

        int iter = OEMSolver::bicgsq(size,op_->systemMatrix(),
            arg.leakPointer(),dest.leakPointer(),epsilon_,maxIter_,verbose_);

        iterations_ = iter;

        std::cout << "OEM-BICGGsq: " << iter << " iterations!\n";
        // finalize operator
        finalize ();
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }

    };


    /** \brief GMRES solver */
    template< class DiscreteFunctionType, class Op >
    class OEMGMRESOp
    : public Operator< DiscreteFunctionType, DiscreteFunctionType >
    {
      typedef OEMGMRESOp< DiscreteFunctionType, Op > This;

    public:
      typedef Op OperatorType;

    private:
      // type of internal projector if no preconditioner given
      typedef OEMSolver :: FakeConditioner FakeConditionerType;

      typedef std::pair < int , double > ReturnValueType;

      template <class OperatorImp, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class DiscreteFunctionImp>
        static ReturnValueType call(OperatorImp & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest,
                         int inner, double eps, int maxIter, bool verbose)
        {
          int size = arg.space().size();
          if(op.hasPreconditionMatrix())
          {
            return OEMSolver::gmres(arg.space().gridPart().comm(),
                     inner,size,op.systemMatrix(),op.preconditionMatrix(),
                     arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
          // in parallel case we need special treatment, if no preconditoner exist
          else if( arg.space().gridPart().comm().size() > 1 )
          {
            OEMSolver::SolverInterfaceImpl< std::decay_t< decltype( op.systemMatrix() ) > > opSolve( op.systemMatrix() );
            FakeConditionerType preConditioner( size, opSolve );
            return OEMSolver::gmres(arg.space().gridPart().comm(),
                     inner,size,op.systemMatrix(),preConditioner,
                     arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
          else
          {
            return OEMSolver::gmres(arg.space().gridPart().comm(),
                     inner,size,op.systemMatrix(),
                     arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
        }
      };

      // without any preconditioning
      template <class OperatorImp>
      struct SolverCaller<OperatorImp,false>
      {
        template <class DiscreteFunctionImp>
        static ReturnValueType call(OperatorImp & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest,
                         int inner, double eps, int maxIter, bool verbose)
        {
          int size = arg.space().size();
          if( arg.space().gridPart().comm().size() > 1 )
          {
            OEMSolver::SolverInterfaceImpl< std::decay_t< decltype( op.systemMatrix() ) > > opSolve( op.systemMatrix() );
            FakeConditionerType preConditioner( size, opSolve );
            return OEMSolver::gmres(arg.space().gridPart().comm(),
                     inner,size,op.systemMatrix(),preConditioner,
                     arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
          else
          {
            return OEMSolver::gmres(arg.space().gridPart().comm(),
                     inner,size,op.systemMatrix(),
                     arg.leakPointer(),dest.leakPointer(),eps,maxIter,verbose);
          }
        }
      };

    public:
      /** \brief constructor of OEM-GMRES
          \param[in] redEps realative tolerance for residual
          \param[in] absLimit absolut solving tolerance for residual
          \param[in] maxIter maximal number of iterations performed
          \param[in] verbose verbosity
      */
      OEMGMRESOp ( double redEps, double absLimit, int maxIter, bool verbose,
                   const ParameterReader &parameter = Parameter::container() )
      : epsilon_( absLimit ),
        maxIter_( maxIter ),
        restart_( parameter.getValue< int >( "oemsolver.gmres.restart", 20 ) ),
        verbose_( verbose ),
        iterations_( 0 )
      {}

      OEMGMRESOp ( double redEps, double absLimit, int maxIter,
                   const ParameterReader &parameter = Parameter::container() )
      : OEMGMRESOp( redEps, absLimit, maxIter,
                    parameter.getValue< bool >( "fem.solver.verbose", false ), parameter )
      {}

      OEMGMRESOp ( double redEps, double absLimit,
                   const ParameterReader &parameter = Parameter::container() )
      : OEMGMRESOp( redEps, absLimit, std::numeric_limits< int >::max(),
                    parameter.getValue< bool >( "fem.solver.verbose", false ), parameter )
      {}

      /** \brief constructor of OEM-GMRES
          \param[in] op Operator to invert
          \param[in] redEps realative tolerance for residual
          \param[in] absLimit absolut solving tolerance for residual
          \param[in] maxIter maximal number of iterations performed
          \param[in] verbose verbosity
      */
      OEMGMRESOp ( OperatorType &op,
                   double redEps, double absLimit, int maxIter, bool verbose,
                   const ParameterReader &parameter = Parameter::container() )
      : OEMGMRESOp( redEps, absLimit, maxIter, verbose, parameter )
      {
        bind( op );
      }

      OEMGMRESOp ( OperatorType &op,
                   double redEps, double absLimit, int maxIter,
                   const ParameterReader &parameter = Parameter::container() )
      : OEMGMRESOp( redEps, absLimit, maxIter, parameter )
      {
        bind( op );
      }

      OEMGMRESOp ( OperatorType &op,
                   double redEps, double absLimit,
                   const ParameterReader &parameter = Parameter::container() )
      : OEMGMRESOp( redEps, absLimit, parameter )
      {
        bind( op );
      }

      void bind ( OperatorType &op )
      {
        op_ = &op;
      }

      void unbind ()
      {
        op_ = nullptr;
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {
      }

      int iterations() const
      {
        return iterations_;
      }

      void setMaxIterations ( int maxIter )
      {
        maxIter_ = maxIter;
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        assert( op_ );
        // prepare operator
        prepare ( arg, dest );

        int size = arg.space().size();
        int inner = std::min( size, restart_ );

        ReturnValueType val =
          SolverCaller<OperatorType,
                       // check wheter operator has precondition methods
                       // to enable preconditioning derive your operator from
                       // OEMSolver::PreconditionInterface
                       std::is_convertible<OperatorType, OEMSolver::PreconditionInterface > ::value >::
                         // call solver, see above
                         call(*op_,arg,dest,inner,epsilon_,maxIter_,verbose_);

        iterations_ = val.first;

        if( arg.space().gridPart().comm().rank() == 0 && verbose_ )
        {
          std::cout << "OEM-GMRES: " << val.first << " iterations! Error: " << val.second << "\n";
        }

        // finalize operator
        finalize ();
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }

    private:
      // no const reference, we make const later
      OperatorType *op_ = nullptr;
      typename DiscreteFunctionType::RangeFieldType epsilon_;
      int maxIter_;
      int restart_;
      bool verbose_ ;
      mutable int iterations_;
    };

    /**
       @}
    **/
#ifdef USE_PARDG_ODE_SOLVER
    /////////////////////////////////////////////////////////////////
    //
    //  GMRES Version of Dennis code
    //
    /////////////////////////////////////////////////////////////////
    // \brief GMRES implementation from Dennis D.
    template <class DiscreteFunctionType, class OperatorImp>
    class GMRESOp : public Operator<
                DiscreteFunctionType,DiscreteFunctionType>
    {
    public:
      typedef OperatorImp OperatorType;
    private:
      typedef OEMSolver :: FakeConditioner FakeConditioner;

      template <class SolverType, bool hasPreconditioning>
      struct SolverCaller
      {
        template< class OperatorImpA, class PreConMatrix, class DiscreteFunctionImp >
        static void solve(SolverType & solver,
                   OperatorImpA & op,
                   const PreConMatrix & pm,
                   const DiscreteFunctionImp & arg,
                   DiscreteFunctionImp & dest)
        {
          int size = arg.space().size();

          OEMSolver::SolverInterfaceImpl<OperatorImpA> opSolve(op,size);

          // in parallel runs we need fake pre conditioner to
          // project vectors onto interior
          if(op.hasPreconditionMatrix())
          {
            OEMSolver::PreconditionerImpl<PreConMatrix> pre(pm,size);
            solver.set_preconditioner(pre);

            // note argument and destination are toggled
            solver.solve(opSolve, dest.leakPointer() , arg.leakPointer() );

            solver.unset_preconditioner();
          }
          else
          {
            // note argument and destination are toggled
            solver.solve(opSolve, dest.leakPointer() , arg.leakPointer() );
          }
        }

        template <class OperatorImpA, class DiscreteFunctionImp>
        static void call(SolverType & solver,
                         OperatorImpA & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest)
        {
          solve(solver,op,op.preconditionMatrix(),arg,dest);
        }
      };

      // without any preconditioning
      template <class SolverType>
      struct SolverCaller<SolverType,false>
      {
        template <class OperatorImpA, class DiscreteFunctionImp>
        static void call(SolverType & solver,
                         OperatorImpA & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest)
        {
          int size = arg.space().size();
          OEMSolver::SolverInterfaceImpl<OperatorImpA> opSolve(op,size);

          // in parallel runs we need fake pre conditioner to
          // project vectors onto interior
          if(arg.space().gridPart().comm().size() > 1)
          {
            FakeConditioner fake(size,opSolve);
            OEMSolver::SolverInterfaceImpl<FakeConditioner> pre(fake);
            solver.set_preconditioner(pre);

            // note argument and destination are toggled
            solver.solve(opSolve, dest.leakPointer() , arg.leakPointer() );
            solver.unset_preconditioner();
          }
          else
          {
            // note argument and destination are toggled
            solver.solve(opSolve, dest.leakPointer() , arg.leakPointer() );
          }
        }
      };

      // solver
      typedef PARDG::GMRES SolverType;
      mutable SolverType solver_;

      // wrapper to fit interface of FGMRES operator
      OperatorType *op_ = nullptr;

      typename DiscreteFunctionType::RangeFieldType epsilon_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;

    public:
      GMRESOp( double  redEps , double absLimit , int maxIter , bool verbose,
          const ParameterReader &parameter = Parameter::container() )
      : solver_(PARDG::Communicator::instance(), parameter.getValue< int >( "oemsolver.gmres.restart", 20 ) )
        , epsilon_ ( absLimit )
        , maxIter_ (maxIter ) , verbose_ ( verbose )
        , iterations_ ( 0 )
      {
      }

      GMRESOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose,
          const ParameterReader &parameter = Parameter::container() )
      : GMRESOp( redEps, absLimit, maxIter, verbose )
      {
        bind( op );
      }

      void bind ( OperatorType &op )
      {
        op_ = &op;
      }

      void unbind ()
      {
        op_ = nullptr;
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {
      }

      int iterations () const
      {
        return iterations_;
      }

      void setMaxIterations ( int maxIter )
      {
        maxIter_ = maxIter;
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        assert( op_ );
        // prepare operator
        prepare ( arg, dest );

        solver_.set_tolerance(epsilon_);
        solver_.set_max_number_of_iterations( maxIter_ );

        if(verbose_)
        {
          solver_.IterativeSolver::set_output(std::cout);
          solver_.DynamicalObject::set_output(std::cout);
        }

        SolverCaller<SolverType,
                       // check wheter operator has precondition methods
                       // to enable preconditioning derive your operator from
                       // OEMSolver::PreconditionInterface
                       std::is_convertible<OperatorType, OEMSolver::PreconditionInterface > ::value >::
                       // call solver, see above
                       call(solver_,op_->systemMatrix(),arg,dest);

        iterations_ = solver_.number_of_iterations();

        // finalize operator
        finalize ();
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }
    };

    template <class DiscreteFunctionType, class OperatorImp>
    class FGMRESOp : public Operator<
                DiscreteFunctionType,DiscreteFunctionType>
    {
    public:
      typedef OperatorImp OperatorType;
    private:
      typedef OEMSolver :: FakeConditioner FakeConditionerType;

      template <class SolverType, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class OperatorImpA, class PreConMatrix, class DiscreteFunctionImp>
        static void solve(SolverType & solver,
                   OperatorImpA & op,
                   const PreConMatrix & pm,
                   const DiscreteFunctionImp & arg,
                   DiscreteFunctionImp & dest)
        {
          int size = arg.space().size();
          OEMSolver::SolverInterfaceImpl<OperatorImpA> opSolve(op,size);
          OEMSolver::PreconditionerImpl<PreConMatrix> pre(pm,size);
          solver.set_preconditioner(pre);

          // note argument and destination are toggled
          solver.solve(opSolve, dest.leakPointer() , arg.leakPointer() );
          solver.unset_preconditioner();
        }

        template <class OperatorImpA, class DiscreteFunctionImp>
        static void call(SolverType & solver,
                         OperatorImpA & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest)
        {
          if(op.hasPreconditionMatrix() )
          {
            solve(solver,op.systemMatrix(),op.preconditionMatrix(),arg,dest);
          }
          else
          {
            SolverCaller<SolverType,false>::call(solver,op,arg,dest);
          }
        }
      };

      // without any preconditioning
      template <class SolverType>
      struct SolverCaller<SolverType,false>
      {
        template <class OperatorImpA, class DiscreteFunctionImp>
        static void solve(SolverType & solver,
                   OperatorImpA & op,
                   const DiscreteFunctionImp & arg,
                   DiscreteFunctionImp & dest)
        {
          int size = arg.space().size();
          OEMSolver::SolverInterfaceImpl<OperatorImpA> opSolve(op,size);
          FakeConditionerType fake(size,opSolve);
          SolverCaller<SolverType,true>::solve(solver,op,fake,arg,dest);
        }

        template <class OperatorImpA, class DiscreteFunctionImp>
        static void call(SolverType & solver,
                         OperatorImpA & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest)
        {
          // not working yet
          assert( false );
          solve(solver,op.systemMatrix(),arg,dest);
        }
      };

      // solver
      typedef PARDG::FGMRES SolverType;
      mutable SolverType solver_;

      // wrapper to fit interface of FGMRES operator
      OperatorType *op_ = nullptr;

      typename DiscreteFunctionType::RangeFieldType epsilon_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;

    public:
      FGMRESOp( double  redEps , double absLimit , int maxIter , bool verbose,
                const ParameterReader &parameter = Parameter::container() )
      : solver_(PARDG::Communicator::instance(), parameter.getValue< int >( "oemsolver.gmres.restart", 20 ) )
        , epsilon_ ( absLimit )
        , maxIter_ (maxIter ) , verbose_ ( verbose )
        , iterations_( 0 )
      {
      }

      FGMRESOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose,
                const ParameterReader &parameter = Parameter::container() )
      : FGMRESOp( redEps, absLimit, maxIter, verbose, parameter )
      {
        bind( op );
      }

      void bind ( OperatorType &op )
      {
        op_ = &op;
      }

      void unbind ()
      {
        op_ = nullptr;
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {
      }

      int iterations () const
      {
        return iterations_;
      }

      void setMaxIterations ( int maxIter )
      {
        maxIter_ = maxIter;
      }

      //! solve the system
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        assert( op_ );
        // prepare operator
        prepare ( arg, dest );

        solver_.set_tolerance(epsilon_);
        solver_.set_max_number_of_iterations( maxIter_ );

        if(verbose_)
        {
          solver_.IterativeSolver::set_output(std::cout);
          solver_.DynamicalObject::set_output(std::cout);
        }

        SolverCaller<SolverType,
                       // check wheter operator has precondition methods
                       // to enable preconditioning derive your operator from
                       // OEMSolver::PreconditionInterface
                       std::is_convertible<OperatorType, OEMSolver::PreconditionInterface > ::value >::
                       // call solver, see above
                       call(solver_,*op_,arg,dest);

        iterations_ = solver_.number_of_iterations();

        // finalize operator
        finalize ();
      }

      //! solve the system
      void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }

    };

    /////////////////////////////////////////////////////////////////
    //
    //  BICGstab Version of Dennis code
    //
    /////////////////////////////////////////////////////////////////
    /*
      \interface
      \brief BICG-stab implementation from Dennis D.
    */
    template <class DiscreteFunctionType, class OperatorImp>
    class BICGSTABOp : public Operator<
                DiscreteFunctionType,DiscreteFunctionType>
    {
    public:
      typedef OperatorImp OperatorType;
    private:
      template <class SolverType, bool hasPreconditioning>
      struct SolverCaller
      {
        template <class OperatorImpA, class PreConMatrix, class DiscreteFunctionImp>
        static void solve(SolverType & solver,
                   OperatorImpA & op,
                   const PreConMatrix & pm,
                   const DiscreteFunctionImp & arg,
                   DiscreteFunctionImp & dest)
        {
          int size = arg.space().size();
          OEMSolver::SolverInterfaceImpl<OperatorImpA> opSolve(op,size);

          OEMSolver::PreconditionerImpl<PreConMatrix> pre(pm,size);
          solver.set_preconditioner(pre);

          // note argument and destination are toggled
          solver.solve(opSolve, dest.leakPointer() , arg.leakPointer() );
          solver.unset_preconditioner();
        }

        template <class OperatorImpA, class DiscreteFunctionImp>
        static void call(SolverType & solver,
                         OperatorImpA & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest)
        {
          if(op.hasPreconditionMatrix())
          {
            solve(solver,op.systemMatrix(),op.preconditionMatrix(),arg,dest);
          }
          else
          {
            SolverCaller<SolverType,false>::call(solver,op,arg,dest);
          }
        }
      };

      // without any preconditioning
      template <class SolverType>
      struct SolverCaller<SolverType,false>
      {
        template <class OperatorImpA, class DiscreteFunctionImp>
        static void solve(SolverType & solver,
                   OperatorImpA & op,
                   const DiscreteFunctionImp & arg,
                   DiscreteFunctionImp & dest)
        {
          int size = arg.space().size();
          OEMSolver::SolverInterfaceImpl<OperatorImpA> opSolve(op,size);

          // note argument and destination are toggled
          solver.solve(opSolve, dest.leakPointer() , arg.leakPointer() );
        }
        template <class OperatorImpA, class DiscreteFunctionImp>
        static void call(SolverType & solver,
                         OperatorImpA & op,
                         const DiscreteFunctionImp & arg,
                         DiscreteFunctionImp & dest)
        {
          solve(solver,op.systemMatrix(),arg,dest);
        }
      };

      // solver
      typedef PARDG::BICGSTAB SolverType;
      mutable SolverType solver_;
      // wrapper to fit interface of GMRES operator
      OperatorType *op_ = nullptr;

      typename DiscreteFunctionType::RangeFieldType epsilon_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;

    public:
      BICGSTABOp( double  redEps , double absLimit , int maxIter , bool verbose ,
                  const ParameterReader &parameter = Parameter::container() )
      : solver_(PARDG::Communicator::instance())
        , epsilon_ ( absLimit )
        , maxIter_ (maxIter ) , verbose_ ( verbose )
        , iterations_( 0 )
      {
      }

      BICGSTABOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose ,
                  const ParameterReader &parameter = Parameter::container() )
      : BICGSTABOp( redEps, absLimit, maxIter, verbose, parameter )
      {
        bind( op );
      }

      void bind ( OperatorType &op )
      {
        op_ = &op;
      }

      void unbind ()
      {
        op_ = nullptr;
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {
      }

      int iterations () const
      {
        return iterations_;
      }

      void setMaxIterations ( int maxIter )
      {
        maxIter_ = maxIter;
      }

      //! solve the system
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        assert( op_ );
        // prepare operator
        prepare ( arg, dest );

        solver_.set_tolerance( epsilon_ );
        solver_.set_max_number_of_iterations( maxIter_ );

        if( verbose_ )
        {
          solver_.IterativeSolver::set_output(std::cout);
          solver_.DynamicalObject::set_output(std::cout);
        }

        SolverCaller<SolverType,
                       // check wheter operator has precondition methods
                       // to enable preconditioning derive your operator from
                       // OEMSolver::PreconditionInterface
                       std::is_convertible<OperatorType, OEMSolver::PreconditionInterface > ::value >::
                       // call solver, see above
                       call(solver_,*op_,arg,dest);

        iterations_ = solver_.number_of_iterations();

        // finalize operator
        finalize ();
      }

      //! solve the system
      void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }

    };
#endif

  } // namespace Fem

} // namespace Dune

#endif //#indef DUNE_FEM_OEMSOLVER_HH
