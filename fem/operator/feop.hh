/**************************************************************************
**       Title: feop
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description: Implementation of a finite element operator class for 
**              general elliptic problems. An example of the use is given
**              in fem/examples/elliptic.
**
**************************************************************************/

#ifndef DUNE_FEOP_HH
#define DUNE_FEOP_HH

//- Dune includes 
#include <dune/common/fmatrix.hh>
#include <dune/grid/common/referenceelements.hh>

//- local includes 
#include "common/operator.hh"
#include "common/localoperator.hh"
#include "feop/spmatrix.hh"
#include "lagrangedofhandler.hh"

namespace Dune {

/*======================================================================*/
/*! @ingroup CGFiniteElement  
 *  \class FEOp 
 *  \brief The FEOp class provides an example of a Finite Element Operator 
 *     
 *  The class is an operator used for general elliptic problems, specialized 
 *  to Lagrange-bases, as explicit boundary setting and kronecker-kills
 *  in the matrix are performed.
 *
 *  The derivation from the FEOpInterface seems superfluous, as the local 
 *  element matrix provider is now specified by a template-parameter. 
 *
 *  The class is used for general elliptic problems + boundary treatment: 
 *  \f{eqnarray*}
 *                 - div(a*grad(u) - b*u) + c*u = f     in Omega    \\
 *                                            u = g_D   in \Gamma_D \\
 *                           (a*grad(u) -b*u) n = g_N   in \Gamma_N \\
 *                 (a*grad(u) -b*u) n + alpha*u = g_R   in \Gamma_R 
 *  \f}
 *
 *  where \l$ a,b,c,g_D,g_N,g_R \l$ are space dependent, alpha a constant and the 
 *  quantities denote
 *              "stiffness"        a(x) 
 *              "velocity"         b(x) 
 *              "mass"             c(x) 
 *              "source"           f(x) 
 *              "dirichletValues"  g_D
 *              "neumannValues"    g_N
 *              "robinValues"      g_R
 *              "alpha"            alpha                   
 *
 * the following assumptions on the basis/model are made:
 *    - The discrete function space is a nodal basis, 
 *      there exists a set of x_i such that phi_j(x_i) = kronecker(i,j)  
 *      the access to these points is done by a LagrangeDOFHandler class.
 *    - a basis function phi_i of a neuman boundary point x_i 
 *      vanishes on the Robin-boundary and vice-versa. A basis function of an 
 *      interior point x_i vanishes on the boundary
 *    - The Dirichlet-Boundary is a closed set, in particular 
 *      point-evaluations indicate whether any of the adjacent edges or 
 *      faces are Dirichlet-boundaries
 *    - complete cell boundaries are of one type (except perhaps their 
 *      lower-codim boundaries): cog-evaluation indicates, whether 
 *      all Lagrange nodes on the intersection are Dirichlet-vertices. 
 *      And cog-evaluation indicates, 
 *      whether whole intersection is of one type, e.g. Neuman or Robin, for 
 *      integration over it. 
 *
 *  weak formulation of the above problem and restriction to the discrete 
 *  function with u_h := sum_j u_j phi_j leads to a system for the 
 *  DOF-vector u:
 *  
 *          M u = b
 *
 *  with
 *               
 *             /   kronecker(i,j)         if x_i is Dirichlet-LagrangePoint 
 *            /
 *    M_ij :=<       \int_\Omega   [a     grad(phi_j) ]^T  grad(phi_i) 
 *            \   -  \int_\Omega   [b     phi_j]^T         grad(phi_i)
 *             \  +  \int_\Omega   c          phi_i       phi_j
 *              \ +  \int_\Gamma_R alpha      phi_i       phi_j      otherwise
 *
 *  and
 *
 *           /    g_D(x_i)               if x_i is Dirichlet-LagrangePoint
 *    b_i :=<   
 *           \      \int_\Omega   f   phi_i
 *            \   + \int_\Gamma_N g_N phi_i
 *             \  + \int_\Gamma_R g_R phi_i                        otherwise
 *
 *  The right hand side is assumed to be assembled by another class, e.g.
 *  RhsAssembler, which is based on element-wise contributions
 *  by a RhsIntegrator class, etc.  
 *
 *  The matrix M has kronecker rows for all dirichlet points, but no 
 *  kronecker-columns. 
 *
 *  Optionally, the matrix and the right hand side can
 *  be processed to have also kronecker-columns, which is beneficial in case
 *  of a symmetric problem. This is done by calling the methods (the latter 
 *  possibly being repeated for changing or multiple rhsides.)
 *
 *       matrixKroneckerColumnsTreatment();
 *       rhsKroneckerColumnsTreatment(rhs);
 * 
 *  This results in
 *
 *                M_sym u = b_sym
 *
 *  The new matrix has entries
 *
 *               /   kronecker(i,j)    if x_i OR x_j is Dirichlet-Lagr.Point 
 *   M_ij_sym :=<    M_ij              otherwise
 * 
 *  The new right hand side has entries:
 *
 *               /    b_i                  if x_i is Dirichlet-Lagr..Point
 *    b_i_sym :=<   
 *               \    b_i - sum_{x_j Dirichlet-Point} M_ij g_D(x_j)   otherwise
 *
 *  In this case, the deleted matrix entries are stored, as they are 
 *  required for every subsequent right hand side modification.
 *
 *  The class depends on two template parameters, a SystemMatrixImp and an 
 *  ElementMatrixIntegratorImp class
 *
 *  The SystemMatrixImp class must provide a storage type for the global 
 *  operator matrix and some arithmetics and basic functionality. In 
 *  particular a print() method and an apply() method are assumed to exist.
 *  Additionally a constructor with three arguments is required as 
 *  explained in allocateSystemMatrix. An add(row,col,val) method is assumed 
 *  to exist. An unitRow(row) and a unitCol(col) method are assumed to exist
 * 
 *  The ElementMatrixIntegratorImp class provides functionality for 
 *  computing a
 *  local element matrix on a given entity without dirichlet-treatment, i.e.
 *  by traversing the grid with the ElementMatrixImp class results in 
 *  the above matrix M_ij, with Dirichlet-values also integrated. After this
 *  Assembly, a Dirichlet-BndCorrection is performed.
 *
 *  Different operating modes are possible (currently however only 
 *  ASSEMBLED is implemented). In case of ASSEMBLED, the whole 
 *  global operator matrix is allocated and completely precomputed by 
 *  corresponding methods. In case of ON_THE_FLY, no complete matrix is
 *  allocated, but the matrix-vector multiplication is performed by 
 *  on-the-fly computation of the elementwise local matrices.
 */
/*======================================================================*/

template <class SystemMatrixImp, class ElementMatrixIntegratorImp>
class FEOp : 
        public Operator<
           typename ElementMatrixIntegratorImp::TraitsType::DiscreteFunctionType::DomainFieldType, 
           typename ElementMatrixIntegratorImp::TraitsType::DiscreteFunctionType::RangeFieldType,
           typename ElementMatrixIntegratorImp::TraitsType::DiscreteFunctionType,
           typename ElementMatrixIntegratorImp::TraitsType::DiscreteFunctionType>

{
  
public:
  //! Type of matrix storage class used for global operator matrix
  typedef SystemMatrixImp SystemMatrixType;  

  //! Type of ElementMatrixContributor
  typedef ElementMatrixIntegratorImp ElementMatrixIntegratorType;  
  
  //! Operation mode: global allocation and multiplication or only 
  //! on-the-fly
  enum OpMode { ON_THE_FLY, ASSEMBLED };

  //! Mode of operation on the matrix/rhs concerning Dirichlet rows/cols
  // enum DirichletTreatmentMode { KRONECKER_ROWS, KRONECKER_ROWS_COLS };

  typedef typename ElementMatrixIntegratorType::TraitsType TraitsType;
  typedef typename TraitsType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename TraitsType::ElementMatrixType ElementMatrixType;
  typedef typename TraitsType::IntersectionQuadratureType 
                   IntersectionQuadratureType;
  typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;

/*======================================================================*/
/*! 
 *   constructor: Initialization of FEOp
 *
 *   Based on an existing instance of an elementmatrixintegrator 
 *   (which knows the model and the function space) the FEOp is initialized. 
 *   Operation mode must be selected as ASSEMBLED or ON_THE_FLY. The number of 
 *   nonzeros in the global matrix must be specified (only relevant in 
 *   ASSEMBLED-Mode). 
 * 
 *   \param elMatInt an instance of an element matrix integrator
 *
 *   \param opMode the operator mode 
 *
 *   \param maxNonZerosPerRow the maximum number of nonzeros per row in the 
 *          global matrix (default 50)
 *
 *   \return the initialized FEOp
 */
/*======================================================================*/
  
  FEOp( ElementMatrixIntegratorType &elMatInt, 
        OpMode opMode = ASSEMBLED,
//        DirichletTreatmentMode dirichletMode = KRONECKER_ROWS,
        int maxNonZerosPerRow = 50) :

          functionSpace_( elMatInt.model().discreteFunctionSpace()),  
          matrix_ (0), 
          matrix_assembled_( false ),
//          arg_ ( NULL ), 
//          dest_ (NULL) , 
          opMode_(opMode),
//          dirichletMode_(dirichletMode),
          maxNonZerosPerRow_(maxNonZerosPerRow),
          elMatInt_(elMatInt),
          isDirichletDOF_(0),
          isDirichletDOF_assembled_(false),
          matrixDirichletColumns_(0)
        {
          // class currently only implemented for non-symmetrized matrix/rhs!
          assert(opMode == ASSEMBLED);
        };
  
/*======================================================================*/
/*! 
 *   destructor: In case of allocation of global operator matrix
 *               it is deallocated, also the dirichletDOF lookup table.
 */
/*======================================================================*/

  ~FEOp( ) 
        {
          if ( matrix_ ) 
              delete matrix_;
          if ( isDirichletDOF_ ) 
              delete isDirichletDOF_;
          if ( matrixDirichletColumns_ ) 
              delete matrixDirichletColumns_;
        };

/*======================================================================*/
/*! 
 *   print: print matrix to standard out, only makes sense in ASSEMBLED mode 
 */
/*======================================================================*/

  void print () const 
  {
    assert(opMode_==ASSEMBLED);
    if(!this->matrix_assembled_) 
        this->assemble();
    this->matrix_->print(std::cout);
  }
 
/*======================================================================*/
/*! 
 *   systemMatrix: return reference to systemMatrix for oem solvers. 
 *
 *   The assembled matrix is returned. That means, if global matrix is not 
 *   yet allocated, a new empty matrix is generated. If current global 
 *   matrix is not yet assembled, an assembly is initiated.
 *   Method only makes sense in ASSEMBLED mode
 *
 *   \return reference to assembled global matrix
 */
/*======================================================================*/
  
  SystemMatrixType& systemMatrix() 
        {
          assert(opMode_==ASSEMBLED);
          //assert(matrix_assembled_ == true);
          if ( !this->matrix_assembled_ )
          {
            if(!this->matrix_)
                allocateSystemMatrix( );
            this->assemble(); 
          }
          return (*this->matrix_);
        }
  
// /*======================================================================*/
// /*! 
//  *   getLocalMatrix: get local element matrix 
//  *
//  *   method providing a local element matrix. This method is required
//  *   for satisfying the FEOpInterface.
//  *
//  *   ??? is this necessary after current implementation of FEOp ???
//  *
//  *   \param 
//  *
//  *   \exception 
//  *
//  *   \return
//  */
// /*======================================================================*/
//   ... getLocalMatrix(...);
  
// /*======================================================================*/
// /*! 
//  *   initialize: delete and deallocate global matrix
//  *
//  *   ??? Why is this method called init, if it does the same as a 
//  *       destructor should do ??? 
//  *
//  */
// /*======================================================================*/
  
//   //! methods for global application of the operator
//   void initialize(){
//     //std::cout << "Matrix reinitialized!" << std::endl ;
    
//     matrix_assembled_ = false;
//     if ( matrix_ ) delete(matrix_);
//     matrix_ = 0;
//   };


/*======================================================================*/
/*! 
 *   operator(): application operator
 *
 *   In case of ASSEMBLED-mode a global
 *   apply() on the matrix is called, i.e. a matrix-vector-multiplication 
 *   with the vector arg is performed, the result stored in dest. This is 
 *   an implicit requirement on the MatrixImp class!
 *   
 *   In case of ON_THE_FLY, the method is not yet implemented.
 *
 *   \param arg reference on the argument for the matrix-vector multiplication 
 *
 *   \param dest reference to storage for the destination vector
 */
/*======================================================================*/

  virtual void operator()(const DiscreteFunctionType &arg, 
                          DiscreteFunctionType &dest) const
  {
    if (opMode_ == ASSEMBLED)
    {
      if (!matrix_assembled_) 
          this->assemble();   
      matrix_->apply( arg, dest );
    }
    else
    {
      std::cout << "operator() in FEOP needs to be implemented for "
                << " ON_THE_FLY! \n";
      assert(1==0);  
    }
  };

// /*======================================================================*/
// /*! 
//  *  prepareGlobal: store argument and destination 
//  *
//  *  The storage for argument and destination is asserted and saved and
//  *  the destination vector is cleared.
//  *
//  *   \param argument and destination storage 
//  */
// /*======================================================================*/
//   void prepareGlobal(const DiscreteFunctionType &Arg, DiscreteFunctionType & Dest )  
//   DUNE_DEPRECATED 
//   { 
//     arg_  = &Arg.argument();
//     dest_ = &Dest.destination();
//     assert(arg_ != NULL); assert(dest_ != NULL);
//     dest_->clear();
//   };

//   //! 

/*======================================================================*/
/*! 
 *   finalize: method required by oemsolvers
 *
 *   method does nothing
 */
/*======================================================================*/
  
  void finalize() 
        { 
// *   The members of the argument and dest-pointer are set to NULL hereby
// *   detaching the user specified storage from any future operations of 
// *   FEOp. 
//     arg_ = NULL; dest_ = NULL;
   };

// ! 
// /*======================================================================*/
// /*! 
//  *   applyLocal: makes local multiply on the fly
//  *
//  *   the global operator matrix is not set up, but the local element 
//  *   matrices are used for incremental matrix-vector multiplication.
//  *   A grid-walkthrough and applyLocal results in a complete matrix-vector
//  *   multiplication
//  * 
//  *   \param the entity
//  */
// /*======================================================================*/

//   template <class EntityType>
//   void applyLocal ( EntityType &en ) const 
//   {
//     const DiscreteFunctionType & arg  = (*arg_);
//     DiscreteFunctionType & dest = (*dest_);

//     typedef typename DiscreteFunctionType::FunctionSpace DiscreteFunctionSpaceType;
//     typedef typename DiscreteFunctionSpaceType::GridType GridType; 
    
//     typedef typename EntityType::IntersectionIterator NeighIt;
    
//     typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType 
//                      BaseFunctionSetType;

//     typedef typename DiscreteFunctionSpaceType::RangeType RangeVecType;
//     typedef typename DiscreteFunctionSpaceType::JacobianRange JacobianRange;
//     typedef typename DiscreteFunctionSpaceType::DomainType DomainVecType;

//     typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
//     typedef typename DiscreteFunctionType::ConstDofIteratorType ConstDofIteratorType;

//     DofIteratorType dest_it = dest.dbegin();
//     ConstDofIteratorType arg_it = arg.dbegin();
      
//     const BaseFunctionSetType & baseSet = functionSpace_.getBaseFunctionSet( en );
//     int numOfBaseFct = baseSet.numBaseFunctions();  

//     // get localMatrixsize from ElementMatrixTraits
//     ...
//     assert( numOfBaseFct <= maxnumOfBaseFct );

//     FieldMatrix<double, maxnumOfBaseFct, maxnumOfBaseFct> mat;
    
//     getLocalMatrix( en, numOfBaseFct, mat);

//     if(this->scalar_ == 1.)
//     {
//       for(int i=0; i<numOfBaseFct; i++) 
//       {  
//         int row = functionSpace_.mapToGlobal( en , i );
//         for (int j=0; j<numOfBaseFct; j++ ) 
//         {
//           int col = functionSpace_.mapToGlobal( en , j );   

//           // scalar comes from LocalOperatorDefault, if operator is scaled,
//           // i.e. with timestepsize
//           dest_it[ row ] += arg_it[ col ] * mat[i][j];
//         }
//       }
//     }
//     else 
//     {
//       for(int i=0; i<numOfBaseFct; i++) 
//       {  
//         int row = functionSpace_.mapToGlobal( en , i );
//         for (int j=0; j<numOfBaseFct; j++ ) 
//         {
//           int col = functionSpace_.mapToGlobal( en , j );   

//           // scalar comes from LocalOperatorDefault, if operator is scaled,
//           // i.e. with timestepsize
//           double val = (this->scalar_) * mat[i][j];

//           dest_it[ row ] += arg_it[ col ] * val;
//         }
//       }
//     }
//   }; // end applyLocal



// /*======================================================================*/
// /*! 
//  *   finalizeLocal: corrects the mapping in order to take into account 
//  *                  dirichlet boundary conditions.
//  * 
//  *   The dofs are assumed to be correlated to the local vertex numbers of 
//  *   the element, so only Lagrange basis are possible by this.
//  *   The dirichlet-treatment is then simple copying of the argument DOF 
//  *   to the destination DOF in case of Dirichlet boundary points. So the 
//  *   assumption is, that the argument already contains correct 
//  *   Dirichlet values.
//  *
//  *   \param the entity on which correction is to be performed
//  */
// /*======================================================================*/
 
//   template <class EntityType>
//   void finalizeLocal ( EntityType &en ) const 
//   {
//     // eliminate the Dirichlet rows and columns 
//     typedef typename DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;
//     typedef typename DiscreteFunctionSpaceType::GridType GridType;
//     typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
//     typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
 
//     const DiscreteFunctionType & arg  = (*arg_);
//     DiscreteFunctionType & dest = (*dest_);

//     const GridPartType & gridPart = arg.getFunctionsSpace().gridPart();

//     typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
//     typedef typename DiscreteFunctionType::ConstDofIteratorType ConstDofIteratorType;

//     DofIteratorType dest_it = dest.dbegin();
//     ConstDofIteratorType arg_it = arg.dbegin();

//     const GeometryType t = en.geometry().type();

//     const IntersectionIteratorType endnit = gridPart.iend(en);
//     for(IntersectionIteratorType nit = gridPart.ibegin(en); nit != endnit ; ++nit)
//     {
//       if(nit.boundary())
//       {
//         int face = nit.numberInSelf();
//         enum { dim = EntityType :: dimension };
//         typedef typename EntityType :: ctype coordType;
        
//         const int faceCodim = 1;

//         if(t.isSimplex())
//         {
//           if( nit.boundaryId() != 0 )
//           {
//             static ReferenceSimplex< coordType, dim > refElem;
//             int novx = refElem.size( face, faceCodim , dim );
//             assert( novx == dim );
//             for(int j=0; j<novx ; j++)
//             {
//               // get all local numbers located on the face 
//               int vx  = refElem.subEntity(face, faceCodim , j , dim );
//               // get global dof numbers of this vertices 
//               int col = functionSpace_.mapToGlobal( en, vx);
//               // set solution on dirichlet bnd 
//               dest_it[col] = arg_it[col];
//             }
//           }
//         }
//         if(t.isCube())
//         {
//           static ReferenceCube< coordType, dim > refElem;
//           int novx = refElem.size( face, faceCodim , dim );
//           for(int j=0; j<novx ; j++)
//           {
//             // get all local numbers located on the face 
//             int vx  = refElem.subEntity(face, faceCodim , j , dim );
//             // get global dof numbers of this vertices 
//             int col = functionSpace_.mapToGlobal( en, vx );
//             // set solution on dirichlet bnd 
//             dest_it[col] = arg_it[col];
//           }
//         }
//       }
//     }
//   }// end finalizeLocal

/*======================================================================*/
/*! 
 *   assemble: perform grid-walkthrough and assemble global matrix
 * 
 *   If the matrix storage is not allocated, new storage is allocated by 
 *   allocateSystemMatrix. 
 *   The begin and end iterators are determined, 
 *   the assembling of the global matrix is initiated by call of 
 *   assembleOnGrid 
 *   and the Dirichlet-rows deleted by bndCorrectMatrixOnGrid. 
 *   The assemled flag is set. 
 *
 *   Method only makes sense in ASSEMBLED-mode
 */
/*======================================================================*/

  void assemble ( ) const
  {
    assert(opMode_==ASSEMBLED);
        
    if(!this->matrix_) 
        allocateSystemMatrix();

    assert(this->matrix_);
    
    {
      // allocate local matrix storage, assumed default constructor on class
      ElementMatrixType mat;      
      
      // generate global matrix without dirichlet-treatment 
      IteratorType it    = functionSpace_.begin(); 
      IteratorType endit = functionSpace_.end(); 

      assembleOnGrid(it, endit, mat);
    }
    
    // generate kronecker-rows for dirichlet DOFs  
    this->bndCorrectMatrix();
    
    // in case of dirichletTreatment by Kronecker-kill: eliminate columns
    //if (dirichletMode_==KRONECKER_ROWS_COLS)
    //    matrixKroneckerColumnsTreatment();
    
    matrix_assembled_ = true;
  };

// the following can be removed, if the RhsAssembly is performed externally

// /*======================================================================*/
// /*! 
//  *   assembleRhs: assembly of rhs of an elliptic problem with dirichlet
//  *                treatment
//  *
//  *   generate the vector b required for solving an elliptic problem M u = b
//  *   the dirichlet-difs are set to g_D.
//  *
//  *   \param rhs a reference to the rhs discrete function
//  */
// /*======================================================================*/
  
//   void assembleRhs(DiscreteFunctionType& rhs)
//         {
//           assert(opMode_==ASSEMBLED);
          
//           typedef typename DiscreteFunctionType::FunctionSpaceType 
//               DiscreteFunctionSpaceType;
//           typedef typename DiscreteFunctionSpaceType::GridType GridType; 
//           typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
          
//           {
//             // generate rhs and generate dirichlet-values in dirichlet DOFs  
//             // both is done in one step, as the BndCorrection can not be done
//             // without grid-walkthrough. So both steps are performed in one
//             // grid walkthrough
//             IteratorType it    = functionSpace_.begin(); 
//             IteratorType endit = functionSpace_.end(); 
            
//             assembleBndCorrectRhsOnGrid(it, endit, rhs);
//           }
          
//           // in case of dirichletTreatment by Kronecker-kill: 
//           // generate vector b_sym;
          
//           if (dirichletMode_==KRONECKER_ROWS_COLS)
//               rhsKroneckerColumnsTreatment(rhs);
//         }
  
/*======================================================================*/
/*! 
 *   markForReassembling: mark the local variables to be no longer actual
 *
 *   method marks the global matrix and the Dirichlet-DOF-list for 
 *   reassembly. This must be called, if the underlying model or 
 *   element-matrix has changed after initialization of the FEOp. By this, 
 *   the next assemble() or apply() or operator() call of the FEOp will 
 *   invoke a new computation of these internal quantities
 */
/*======================================================================*/

  void markForReassembling()
        {
          isDirichletDOF_assembled_ = false;
          matrix_assembled_ = false;
        };

/*======================================================================*/
/*! 
 *   matrixKroneckerColumnsTreatment: produce kronecker-columns in matrix
 *
 *   This method can be called after Matrix assembly.
 *   The method changes the matrix M (given above), such that 
 *   kronecker-columns 
 *   are produced for dirichlet DOFs. For symmetric problems this results in a 
 *   symmetric matrix, which is beneficial for system-solvers. Note, that the 
 *   right-hand side similarly has to be modified by 
 *   rhsKroneckerColumnsTreatment(rhs) in 
 *   order to give the identical result. The old matrix values are stored in 
 *   a member variable in order to be able to process arbitrary many RHS 
 *   vectors afterwards.
 *
 *   The new matrix has entries
 *
 *               /   kronecker(i,j)     if x_i OR x_j is Dirichlet-Bnd.Point 
 *   M_ij_sym :=<    M_ij               otherwise
 */
/*======================================================================*/

  void matrixKroneckerColumnsTreatment() const
        {
          assert(isDirichletDOF_);
          assert(matrix_);
          
          // delete storage if allocated
          if (matrixDirichletColumns_)
              delete(matrixDirichletColumns_);
          
          // count new number of Dirichlet-DOFs
          int nDirDOFs = 0;
          for (int i=0; i!=matrix_->rows(); i++)
              if ((*isDirichletDOF_)(0,i)) // so is Dirichlet-DOF
                  nDirDOFs++;
          
          // allocate new storage for to be deleted matrix entries
          matrixDirichletColumns_ = new SparseRowMatrix<double>
              (matrix_->rows(), 
               matrix_->cols(), 
               nDirDOFs);
          
          assert(matrixDirichletColumns_);
          
          // save matrix entries 
          matrixDirichletColumns_->clear();
          for (int i=0; i!=matrix_->rows(); i++)
          {
            SparseRowMatrix<int>::ColumnIterator it = 
                isDirichletDOF_->rbegin(0);
            for (;it!=isDirichletDOF_->rend(0);++it)
                matrixDirichletColumns_->set(i,it.col(),
                                            (*matrix_)(i,it.col()));
          }
          
          // delete matrix columns
          SparseRowMatrix<int>::ColumnIterator it = 
              isDirichletDOF_->rbegin(0);
          for (;it!=isDirichletDOF_->rend(0);++it)
              matrix_->unitCol(it.col());          
        };

/*======================================================================*/
/*! 
 *   rhsKroneckerColumnsTreatment: modify computed right-hand side for 
 *                               dirichlet nodes.
 *
 *   Method is an optional postprocessing step 
 *   method changes the right-hand-side b (assumed to be given as defined 
 *   above), such that the symmetrized matrix (after application of 
 *   matrixKroneckerColumnsTreatment) 
 *   with the new rhs produces the same solution as the old matrix M with 
 *   the old rhs b. The new vector has entries:
 *
 *               /    b_i                         if x_i is Dirichlet-Bnd.Point
 *    b_i_sym :=<   
 *               \    b_i - sum_{x_j Dirichlet-Point} M_ij g_D(x_j)   otherwise
 *
 *   Note that the values g_D(x_j) are exactly assumed to be the values of
 *   the provided vector b_j, which is used for implementation instead of
 *   reevaluating the boundary values. 
 *
 *   Current implementation is quite slow, should be improved by
 *   rewriting as matrix-vector multiplication.
 *
 *   \param rhs a reference to the assembled rhs b
 */
/*======================================================================*/

  void rhsKroneckerColumnsTreatment(DiscreteFunctionType& rhs) const
        {
          assert(isDirichletDOF_);
          assert(matrixDirichletColumns_);
          assert ( rhs.size() == isDirichletDOF_->cols() );
          
          std::cout << "entered rhsKroneckerColumnsTreatment\n";

          // temporary store current values rhs[i] for Dirichlet-DOFs
          SparseRowMatrix<double> 
              rhsDirichletValues(1,isDirichletDOF_->cols(),
                                 isDirichletDOF_->numNonZeros());
          rhsDirichletValues.clear();          
          DofIteratorType dit = rhs.dbegin();
          for (int i=0; i!=matrix_->rows(); i++, ++dit)
              if ((*isDirichletDOF_)(0,i)) // so is Dirichlet-DOF
                  rhsDirichletValues.set(0,i,*dit);
//          std::cout << "temporarily stored current-Rhs values\n";
          
          // modify values of rhs-vector for non-dirichlet-DOFs 
          dit = rhs.dbegin();
          for (int i=0; i!=matrix_->rows(); i++, ++dit)
              if ((*isDirichletDOF_)(0,i)==0.0) // so is non-Dirichlet
              {
                std::cout << "found non-Dirichlet-DOF " << i 
                          <<" for modification\n";
                SparseRowMatrix<int>::ColumnIterator it = 
                    isDirichletDOF_->rbegin(0);
                const SparseRowMatrix<int>::ColumnIterator endit = 
                    isDirichletDOF_->rend(0);
//                std::cout << "initialized column iterator \n";         
                for (;it!=endit;++it) 
                {
//                  std::cout << "entered DOF-setting loop with " 
//                            << "column iterator  it pointing to col = " 
//                            << it.col() 
//                            << "   \n";         
                    *dit -= 
                        (*matrixDirichletColumns_)(i,it.col())
                        * rhsDirichletValues(0,it.col()); 
//                  std::cout << "setting of  Dirichlet DOF finished\n";         
                }
              }
//          std::cout << "finished Rhs-KroneckerColumnTreatment\n";
        }
  
// private methods used by the public ones.
private:
    
/*======================================================================*/
/*! 
 *   allocateSystemMatrix: allocation of a new global matrix
 *
 *   The SystemMatrixclass is required to have a constructor with the syntax
 *   SystemMatrixType(nrows, ncols, nonzeros_per_row). Deallocation of the 
 *   global matrix is performed in the destructor. 
 *   
 *   Currently the nonzeros per row must be specified in the constructor of FEOp. 
 *   This could be improved, e.g. by a Traitsclass as template-parameter, which 
 *   contains the SystemMatrixType and the maxnumbernonzeros
 */
/*======================================================================*/

  void allocateSystemMatrix( ) const  
        { 
          
          // the following choice of numnonzeros seems strange: must be exponential 
          // with dimension
          // typedef typename DiscreteFunctionType::FunctionSpaceType::GridType GridType; 
          // enum { dim = GridType::dimension };
          // int maxNonZerosPerRow_ = 15 * (dim-1);
          
          // SystemMatrixType* 
          matrix_ =
              new SystemMatrixType( 
                  this->functionSpace_.size ( ) , 
                  this->functionSpace_.size ( ) , 
                  maxNonZerosPerRow_);
          assert(matrix_);
        };
    
/*======================================================================*/
/*! 
 *   assembleOnGrid: perform grid walkthrough and assemble matrix
 *
 *   For each element, the local element matrix is determined into the 
 *   given local matrix storage and distributed into the global matrix.
 *   Distribution is performed by an add(row,col,val) method on the 
 *   global matrix class.
 *
 *   \param start and end iterator and storage for a local matrix 
 */
/*======================================================================*/

  template <class GridIteratorType, class ElementMatrixImp>
  void assembleOnGrid ( GridIteratorType &it, GridIteratorType &endit, 
                              ElementMatrixImp &mat) const
        {
          typedef typename 
              DiscreteFunctionType::FunctionSpaceType::BaseFunctionSetType 
              BaseFunctionSetType;
          
          // run through grid and add up local contributions
          for( ; it != endit; ++it )
          {
            const BaseFunctionSetType & baseSet 
                = this->functionSpace_.getBaseFunctionSet( *it );
            const int numOfBaseFct = baseSet.numBaseFunctions();  
            
            // setup local element matrix 
            // (size check is performed in elementmatrix construction)
            mat.clear();            
            elMatInt_.addElementMatrix( *it, mat, 1.0);
            
            for(int i=0; i<numOfBaseFct; i++) 
            { 
              int row = functionSpace_.mapToGlobal( *it , i );
              for (int j=0; j<numOfBaseFct; j++ ) 
              {
                int col = functionSpace_.mapToGlobal( *it , j );    
                matrix_->add( row , col , mat(i,j));
              }
            }
          }
        }

/*======================================================================*/
/*! 
 *   searchDirichletDOFs: determine lookup table for dirichlet-values
 *
 *   method fills the local vector isDirichletDOF_ as multiple operations with
 *   Dirichlet-boundaries are necessary, e.g. matrix-boundary treatment, 
 *   NOT rhs-assembly but later symmetrization of the system. By this multiple 
 *   global grid walkthroughs can be prevented and replaced by single run 
 *   over the lookup-table.
 */
/*======================================================================*/
  
  void searchDirichletDOFs() const
        {
          typedef typename 
              DiscreteFunctionType::FunctionSpaceType DiscreteFunctionSpaceType;
          typedef typename DiscreteFunctionSpaceType::GridType GridType; 
          typedef typename GridType::template Codim<0>::Entity EntityType;
          typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
          typedef typename GridPartType::IntersectionIteratorType 
              IntersectionIteratorType;
          typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;

          // allocate isDirichlet-Vector if not already done:
          // to be sure: worst case: all DOFS are Dirichlet... undoubtedly 
          // memory waste in highly resolved case, but no idea... 
          // give maximum in Model? Or where?
          if (!this->isDirichletDOF_)
              isDirichletDOF_ = new SparseRowMatrix<int>
                  (1, 
                   functionSpace_.size(),
                   functionSpace_.size());
          
          assert(this->isDirichletDOF_);
          
          // start by filling all DOFs with zero
          isDirichletDOF_->clear(); 
          
          IteratorType it    = functionSpace_.begin(); 
          IteratorType endit = functionSpace_.end(); 
          
          GridPartType & gridPart = functionSpace_.gridPart();
          
          for( ; it != endit; ++it ) 
          {
            const EntityType & en = *it; 
            
            const GeometryType t = en.geometry().type();
            const IntersectionIteratorType endnit = gridPart.iend(en);
            for(IntersectionIteratorType nit = gridPart.ibegin(en); 
                nit != endnit ; ++nit)
            {
              if(nit.boundary())
              {
                // if boundary, get cog of intersection and check for dirichlet
                IntersectionQuadratureType 
                    iquad(nit,0,
                          TraitsType::IntersectionQuadratureType::INSIDE);
                // check, whether real cog integration is performed
                assert(iquad.nop()==1);
                if (elMatInt_.model().boundaryType(en,iquad,0) == 
                    TraitsType::Dirichlet)
                {
                  // determine all local DOF numbers of intersection vertices

//                  typedef typename EntityType :: ctype coordType;
//                  enum { dim = EntityType :: dimension };

//                  ReferenceElementContainer<coordType,dim> refElemCont;
//                  const ReferenceElement<coordType,dim>& 
//                      refElem = refElemCont(t); // t is geometrytype
                  
                  LagrangeDofHandler<DiscreteFunctionSpaceType> 
                      dofHandler(functionSpace_,en);
                  
                  const int faceCodim = 1;
                  int face = nit.numberInSelf();
//                  int novx = refElem.size( face, faceCodim , dim );
                  int novx = dofHandler.numDofsOnFace( face, faceCodim );
                  //  assert( novx == dim );

                  for(int j=0; j<novx ; j++)
                  {
                    // get all local numbers located on the face 
                    // int vx  = refElem.subEntity(face, faceCodim , j , dim );
                    int vx  = dofHandler.entityDofNum(face, faceCodim , j );

                    // get global dof numbers of this vertices 
                    int row = functionSpace_.mapToGlobal( en, vx);
                    // store DOF as DirichletDOF
                    isDirichletDOF_->set(0,row,1);                    
                  }
                }
              }
            }
          }

          isDirichletDOF_assembled_ = true;
          
        }; // end of searchDirichletDOFs 
  
/*======================================================================*/
/*! 
 *   bndCorrectMatrix: treatment of Dirichlet-DOFS
 *
 *   delete rows for dirichlet DOFS, setting diagonal 
 *   element to 1. This is reasonable, as Lagrange Basis is implicitly 
 *   assumed, 
 *   so the RHS being the exact dirichlet-values, therefore, the matrix 
 *   row must be a unit-row.
 */
/*======================================================================*/
  
//  template <class GridIteratorType>
  void bndCorrectMatrix() const
        {
          if (!isDirichletDOF_assembled_)
              searchDirichletDOFs();
          
          // eliminate the Dirichlet rows by converting to unit-rows      
          typename SparseRowMatrix<int>::ColumnIterator 
              it = isDirichletDOF_->rbegin(0); 
          typename SparseRowMatrix<int>::ColumnIterator 
              endit = isDirichletDOF_->rend(0);
          
          for (;it!=endit;++it)
              matrix_->unitRow(it.col());          
        }; // end of bndCorrectMatrixOnGrid

// /*======================================================================*/
// /*! 
//  *   assembleBndCorrectRhsOnGrid: assemble righ hand side vector 
//  *
//  *   the right hand side vector b defined by the current model is computed
//  *   The dirichlet-DOFs are set to exact dirichlet-values.
//  *
//  *   \param it begin iterator for grid walkthrough
//  *
//  *   \param endit end iterator for grid walkthrough
//  *
//  *   \param rhs right hand side vector storage
//  */
// /*======================================================================*/

//   void assembleBndCorrectRhsOnGrid(it, endit, rhs)
//         {
//           ...
//         };
  
//! member variables:
private: 
  
  //! the corresponding function_space 
  DiscreteFunctionSpaceType & functionSpace_;
   
  //! The following member variables will be modified by "operator() const" so 
  //! must be mutable

  //! pointer to the representing global matrix 
  mutable SystemMatrixType *matrix_ ; 

  //! flag indicating whether the global matrix is assembled 
  mutable bool matrix_assembled_;

  //! vector (matrix with single row), which indicates the Dirichlet-DOFs 
  //! This type is made explicit, as not much optimization can be performed 
  //! here, or? 
  mutable SparseRowMatrix<int> *isDirichletDOF_;

  //! flag indicating whether the DirichletDOF lookup table is assembled 
  mutable bool isDirichletDOF_assembled_;

  //! matrix for storage of the matrix columns, which are deleted 
  //! during symmetrization, but required for RHS modification 
  //! This type is made explicit, as not much optimization can be performed 
  //! here, or? 
  mutable SparseRowMatrix<double> *matrixDirichletColumns_;

  //! reference to element matrix generator provided during initialization
  ElementMatrixIntegratorType& elMatInt_;
   
  //! operator mode 
  OpMode opMode_;

  //! dirichlet-treatment mode 
  // DirichletTreatmentMode dirichletMode_;

  //! maximal number of nonzeros per row in the global matrix. Used for 
  //! allocation. 
  int  maxNonZerosPerRow_;

  // //! pointers to storage of argument and destination, only required in 
  // LocalOperator Mode
  // const DiscreteFunctionType * arg_;
  // DiscreteFunctionType * dest_;

};

} // end namespace


#endif
