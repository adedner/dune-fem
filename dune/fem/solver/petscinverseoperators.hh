#ifndef DUNE_FEM_PETSCINVERSEOPERATORS_HH
#define DUNE_FEM_PETSCINVERSEOPERATORS_HH

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/inverseoperatorinterface.hh>

#if HAVE_PETSC
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/solver/parameter.hh>

namespace Dune
{

  namespace Fem
  {

    //=====================================================================
    // Implementation for PETSc matrix based Krylov solvers
    //=====================================================================

    /** @ingroup OEMSolver
        @{
    **/

    // PETScSolver
    // --------------
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
    class PetscInverseOperator;

    template <class DF, class Op >
    struct PetscInverseOperatorTraits
    {
    private:
      typedef typename DF :: DiscreteFunctionSpaceType SpaceType ;
    public:
      typedef DF                                   DiscreteFunctionType;
      typedef Op                                   OperatorType;
      typedef OperatorType                         PreconditionerType;
      typedef PetscDiscreteFunction< SpaceType >   SolverDiscreteFunctionType;
      typedef PetscLinearOperator< DF, DF >        AssembledOperatorType;
      typedef PetscInverseOperator< DF, Op >       InverseOperatorType;
    };

    /** \brief PETSc KSP solver context for PETSc Mat and PETSc Vec */
    template< class DF, class Op >
    class PetscInverseOperator : public InverseOperatorInterface< PetscInverseOperatorTraits< DF, Op > >
    {
    protected:
      // monitor function for PETSc solvers
      static PetscErrorCode
      monitor (KSP ksp, PetscInt it, PetscReal rnorm, void *mctx)
      {
        if( Parameter :: verbose () )
        {
          std::cout << "PETSc::KSP:  it = " << it << "    res = " << rnorm << std::endl;
        }
        return PetscErrorCode(0);
      }

      // destroy solver context
      struct KSPDeleter
      {
        void operator() ( KSP* p ) const
        {
          if( !p )
            return;

          ::Dune::Petsc::KSPDestroy( p );
          delete p;
        }
      };

      typedef PetscInverseOperatorTraits< DF, Op > Traits;
      typedef InverseOperatorInterface< Traits > BaseType;
      friend class InverseOperatorInterface< Traits >;
    public:

      typedef typename BaseType :: SolverDiscreteFunctionType    PetscDiscreteFunctionType;
      typedef typename BaseType :: OperatorType                  OperatorType;

      /** \brief constructor
       *
       *  \param[in] op Mapping describing operator to invert
       *  \param[in] reduction reduction epsilon
       *  \param[in] absLimit absolute limit of residual (not used here)
       *  \param[in] maxIter maximal iteration steps
       *  \param[in] verbose verbosity
       *
       *  \note PETSc Krylov solvers uses the relative reduction.
       */
      [[deprecated]]
      PetscInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             int maxIter,
                             bool verbose,
                             const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {
        parameter_->setLinAbsTol( absLimit );
        parameter_->setLinReduction( reduction );
        parameter_->setMaxLinearIterations( maxIter );
        parameter_->setVerbose( verbose );
        bind( op );
      }

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolute limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      [[deprecated]]
      PetscInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             int maxIter,
                             const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : PetscInverseOperator( reduction, absLimit, maxIter, parameter.verbose(), parameter )
      {
        parameter_->setLinAbsTol( absLimit );
        parameter_->setLinReduction( reduction );
        parameter_->setMaxLinearIterations( maxIter );
        bind( op );
      }

      [[deprecated]]
      PetscInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {
        parameter_->setLinAbsTol( absLimit );
        parameter_->setLinReduction( reduction );
        bind( op );
      }

      /** \brief constructor
       *
       *  \param[in] op Mapping describing operator to invert
       *  \param[in] reduction reduction epsilon
       *  \param[in] absLimit absolute limit of residual (not used here)
       *  \param[in] maxIter maximal iteration steps
       *  \param[in] verbose verbosity
       *
       *  \note PETSc Krylov solvers uses the relative reduction.
       */
      [[deprecated]]
      PetscInverseOperator ( double reduction, double absLimit, int maxIter, bool verbose,
                             const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {
        parameter_->setLinAbsTol( absLimit );
        parameter_->setLinReduction( reduction );
        parameter_->setMaxLinearIterations( maxIter );
        parameter_->setVerbose( verbose );
      }

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolute limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      [[deprecated]]
      PetscInverseOperator ( double reduction, double absLimit, int maxIter,
                             const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {
        parameter_->setLinReduction( reduction );
        parameter_->setLinAbsTol( absLimit );
        parameter_->setMaxLinearIterations( maxIter );
      }

      [[deprecated]]
      PetscInverseOperator ( double reduction, double absLimit,
                             const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {
        parameter_->setLinReduction( reduction );
        parameter_->setLinAbsTol( absLimit );
      }

      [[deprecated]]
      PetscInverseOperator ( double reduction, double absLimit,
                             unsigned int maxIter, bool verbose,
                             const ParameterReader& parameter )
      : BaseType( SolverParameter( parameter ) )
      {
        parameter_->setLinAbsTol( absLimit );
        parameter_->setLinReduction( reduction );
        parameter_->setMaxLinearIterations( maxIter );
        parameter_->setVerbose( verbose );
      }

      //non-deprecated constructors
      PetscInverseOperator ( const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {}

      PetscInverseOperator (  const OperatorType &op, const SolverParameter &parameter = SolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {
        bind( op );
      }

      void bind ( const OperatorType &op )
      {
        BaseType :: bind( op );
        initialize( *parameter_ );
      }

      void unbind ()
      {
        BaseType :: unbind();
        ksp_.reset();
      }

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: " << solverName_ << " eps = " << parameter_->linReduction() ;
        out  << "\\\\ \n";
      }

    protected:
      void initialize ( const SolverParameter& parameter )
      {
        if( !assembledOperator_ )
          DUNE_THROW(NotImplemented,"Petsc solver with matrix free implementations not yet supported!");

        // Create linear solver context
        ksp_.reset( new KSP() );
        ::Dune::Petsc::KSPCreate( &ksp() );

        // attach Matrix to linear solver context
        Mat& A = const_cast<Mat &> (assembledOperator_->petscMatrix());
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5
        ::Dune::Petsc::KSPSetOperators( ksp(), A, A, SAME_PRECONDITIONER);
#else
        ::Dune::Petsc::KSPSetOperators( ksp(), A, A );
#endif

        // allow for non-zero initial guess
        ::Dune::Petsc::KSPSetInitialGuessNonzero( ksp(), PETSC_TRUE );

        // set prescribed tolerances
        PetscInt  maxits = parameter_->maxLinearIterations();
        PetscReal reduc  = parameter_->linReduction();
        ::Dune::Petsc::KSPSetTolerances(ksp(), reduc, 1.e-50, PETSC_DEFAULT, maxits);

        enum class PetscSolver {
            cg        = SolverParameter::cg,
            bicgstab  = SolverParameter::bicgstab,
            gmres     = SolverParameter::gmres,
            minres    = SolverParameter::minres,
            bicg      = minres+1,
            only_prec = bicg+1,
            defaults  = only_prec+1
          };

        // if special petsc solver parameter exists use that one, otherwise
        // use krylovMethod from SolverParameter
        const auto& reader = parameter.parameter();
        PetscSolver kspType = PetscSolver::defaults;
        const std::string kspNames[] = { "cg" , "bicgstab", "gmres", "minres", "bicg", "preonly", "defaults" };
        if( reader.exists("petsc.kspsolver.method") )
        {
          // see PETSc docu for more types
          kspType = static_cast< PetscSolver >( reader.getEnum("petsc.kspsolver.method", kspNames, int(PetscSolver::defaults) ) );
        }
        else
          kspType = static_cast< PetscSolver >( parameter.krylovMethod() );

        solverName_ = kspNames[ static_cast< int >( kspType ) ];

        //  select linear solver
        switch( kspType )
        {
          case PetscSolver::cg:
            ::Dune::Petsc::KSPSetType( ksp(), KSPCG );
            break;
          case PetscSolver::bicgstab:
            ::Dune::Petsc::KSPSetType( ksp(), KSPBCGS );
              break;
          case PetscSolver::gmres:
            {
              ::Dune::Petsc::KSPSetType( ksp(), KSPGMRES );
              PetscInt restart = 10;
              if( reader.exists("petsc.gmresrestart") )
              {
                restart = reader.getValue<int>("petsc.gmresrestart", restart );
              }
              else
                restart = parameter.gmresRestart() ;

              ::Dune::Petsc::KSPGMRESSetRestart( ksp(), restart );
              break;
            }
          case PetscSolver::minres:
            ::Dune::Petsc::KSPSetType( ksp(), KSPMINRES );
            break;
          case PetscSolver::bicg:
            ::Dune::Petsc::KSPSetType( ksp(), KSPBICG );
              break;
          case PetscSolver::only_prec:
            ::Dune::Petsc::KSPSetType( ksp(), KSPPREONLY );
              break;
          case PetscSolver::defaults:
            // setup solver context from database/cmdline options
            ::Dune::Petsc::KSPSetFromOptions( ksp() );
            ::Dune::Petsc::KSPSetUp( ksp() );
            break;
          default:
            DUNE_THROW(InvalidStateException,"PetscInverseOperator: invalid solver choosen." );
        }

        /////////////////////////////////////////////
        //  preconditioning
        /////////////////////////////////////////////

        enum class PetscPrec {
            defaults  = 0,
            none      = 1,   // no preconditioning
             // parallel preconditioners
            oas       = 2,   // Overlapping Additive Schwarz
            sor       = 3,   // SOR and SSOR
            jacobi    = 4,   // Jacobi preconditioning
             // requiring additional packages
            hypre     = 5,   // Hypre preconditioning
            ml        = 6,   // ML preconditioner (from Trilinos)
            // serial preconditioners
            ilu       = 7,   // ILU preconditioning
            icc       = 8,   // Incomplete Cholesky factorization
            // direct solvers
            lu        = 9,   // LU factorization
          };

        const std::string pcNames[] = { "default", "none", "asm", "sor", "jacobi", "hypre", "ml", "ilu", "icc", "lu" };
        PetscPrec pcType = static_cast< PetscPrec >( reader.getEnum("petsc.preconditioning.method", pcNames, int(PetscPrec::defaults) ) );

        // setup preconditioning context
        PC pc;
        ::Dune::Petsc::KSPGetPC( ksp(), &pc );

        switch( pcType )
        {
          case PetscPrec::defaults:
            // don't setup the pc context twice
            if ( kspType != PetscSolver::defaults )
            {
              // setup pc context from database/cmdline options
              ::Dune::Petsc::PCSetFromOptions( pc );
              ::Dune::Petsc::PCSetUp( pc );
            }
            break;
          case PetscPrec::none:
            ::Dune::Petsc::PCSetType( pc, PCNONE );
            break;
          case PetscPrec::oas:
            {
              ::Dune::Petsc::PCSetType( pc, PCASM );
              ::Dune::Petsc::PCSetUp( pc );
              break;
            }
          case PetscPrec::sor:
            ::Dune::Petsc::PCSetType( pc, PCSOR );
            break;
          case PetscPrec::jacobi:
            ::Dune::Petsc::PCSetType( pc, PCJACOBI );
            break;
          case PetscPrec::hypre:
            {
#ifdef PETSC_HAVE_HYPRE
              enum class HyprePrec {
                  boomeramg = 0,
                  parasails = 1,
                  pilut = 2
                };

              const std::string hyprePCNames[] = { "boomer-amg", "parasails", "pilu-t" };
              HyprePrec hypreType = static_cast< HyprePrec >( reader.getEnum( "petsc.preconditioning.hypre.method", hyprePCNames, 0 ) );

              std::string hypre;
              if ( hypreType == HyprePrec::boomeramg )
                hypre = "boomeramg";
              else if ( hypreType == HyprePrec::parasails )
                hypre = "parasails";
              else if ( hypreType == HyprePrec::pilut )
                hypre = "pilut";
              else
                DUNE_THROW( InvalidStateException, "PetscInverseOperator: invalid hypre preconditioner choosen." );

              ::Dune::Petsc::PCSetType( pc, PCHYPRE );
              ::Dune::Petsc::PCHYPRESetType( pc, hypre.c_str() );
              ::Dune::Petsc::PCSetUp( pc );
#else // PETSC_HAVE_HYPRE
              DUNE_THROW( InvalidStateException, "PetscInverseOperator: petsc not build with hyper support." );
#endif // PETSC_HAVE_HYPRE
              break;
            }
          case PetscPrec::ml:
#ifdef PETSC_HAVE_ML
            ::Dune::Petsc::PCSetType( pc, PCML );
#else // PETSC_HAVE_ML
              DUNE_THROW( InvalidStateException, "PetscInverseOperator: petsc not build with ml support." );
#endif // PETSC_HAVE_ML
            break;
          case PetscPrec::ilu:
            {
              if ( MPIManager::size() > 1 )
                DUNE_THROW( InvalidStateException, "PetscInverseOperator: ilu preconditioner does not work in parallel." );

              // get fill-in level
              PetscInt pcLevel = reader.getValue<int>("petsc.preconditioning.levels", 0 );

              ::Dune::Petsc::PCSetType( pc, PCILU );
              ::Dune::Petsc::PCFactorSetLevels( pc, pcLevel );
              break;
            }
            ::Dune::Petsc::PCSetType( pc, PCML );
          case PetscPrec::icc:
            {
#ifdef PETSC_HAVE_ICC
              if ( MPIManager::size() > 1 )
                DUNE_THROW( InvalidStateException, "PetscInverseOperator: icc preconditioner does not worl in parallel." );

              // get fill-in level
              PetscInt pcLevel = reader.getValue<int>("petsc.preconditioning.levels", 0 );

              ::Dune::Petsc::PCSetType( pc, PCICC );
              ::Dune::Petsc::PCFactorSetLevels( pc, pcLevel );
#else // PETSC_HAVE_ICC
              DUNE_THROW( InvalidStateException, "PetscInverseOperator: petsc not build with icc support." );
#endif // PETSC_HAVE_ICC
              break;
            }
          case PetscPrec::lu:
            {
              enum class Factorization {
                  petsc = 0,
                  superlu = 1,
                  mumps = 2
                };

              const std::string factorizationNames[] = { "petsc", "superlu", "mumps" };
              Factorization factorType = static_cast< Factorization >( reader.getEnum( "petsc.preconditioning.lu.method", factorizationNames, 0 ) );

              ::Dune::Petsc::PCSetType( pc, PCLU );

              if ( factorType == Factorization::petsc )
                ::Dune::Petsc::PCFactorSetMatSolverPackage( pc, MATSOLVERPETSC );
              else if ( factorType == Factorization::superlu )
                ::Dune::Petsc::PCFactorSetMatSolverPackage( pc, MATSOLVERSUPERLU_DIST );
              else if ( factorType == Factorization::mumps )
                ::Dune::Petsc::PCFactorSetMatSolverPackage( pc, MATSOLVERMUMPS );
              else
                DUNE_THROW( InvalidStateException, "PetscInverseOperator: invalid factorization package choosen." );

              ::Dune::Petsc::PCSetUp( pc );
              break;
            }
          default:
            DUNE_THROW( InvalidStateException, "PetscInverseOperator: invalid preconditioner choosen." );
        }

        // set monitor in verbose mode
        if( parameter.verbose() )
        {
          ::Dune::Petsc::KSPView( ksp() );
          ::Dune::Petsc::KSPMonitorSet( ksp(), &monitor, PETSC_NULL, PETSC_NULL);
        }
      }

      int apply( const PetscDiscreteFunctionType& arg, PetscDiscreteFunctionType& dest ) const
      {
        // need to have a 'distributed' destination vector for continuous spaces
        if( dest.space().continuous() )
          dest.dofVector().clearGhost();

        // call PETSc solvers
        ::Dune::Petsc::KSPSolve( *ksp_, *arg.petscVec() , *dest.petscVec() );

        // a continuous solution is 'distributed' so need a communication here
        if( dest.space().continuous() )
        {
          dest.communicate();
        }

        // get number of iterations
        PetscInt its ;
        ::Dune::Petsc::KSPGetIterationNumber( *ksp_, &its );
        KSPConvergedReason reason;
        ::Dune::Petsc::KSPGetConvergedReason( *ksp_, &reason );

        bool converged = int(reason) >= 0 ;

        return (converged) ? its : -its;
      }

    protected:
      KSP & ksp () { assert( ksp_ ); return *ksp_; }

      using BaseType :: assembledOperator_;
      using BaseType :: parameter_;
      using BaseType :: iterations_;

      std::unique_ptr< KSP, KSPDeleter > ksp_;   // PETSc Krylov Space solver context

      std::string solverName_;
    };

  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCINVERSEOPERATORS_HH
