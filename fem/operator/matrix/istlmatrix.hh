#ifndef DUNE_ISTLMATRIXWRAPPER_HH
#define DUNE_ISTLMATRIXWRAPPER_HH

#if HAVE_DUNE_ISTL

//- system includes 
#include <vector> 

//- Dune istl includes 
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>

//- Dune fem includes 
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/operator/common/localmatrix.hh>

namespace Dune { 

  ///////////////////////////////////////////////////////
  // --BlockMatrixHandle
  //////////////////////////////////////////////////////
  template <class LittleBlockType, class DiscreteFunctionType> 
  class ImprovedBCRSMatrix : public BCRSMatrix<LittleBlockType> 
  {
    public:
      typedef BCRSMatrix<LittleBlockType> BaseType; 
      typedef typename BaseType :: RowIterator RowIteratorType ;
      typedef typename BaseType :: ColIterator ColIteratorType ;

      typedef typename BaseType :: size_type size_type;

      //===== type definitions and constants

      //! export the type representing the field
      typedef typename BaseType::field_type field_type;

      //! export the type representing the components
      typedef typename BaseType::block_type block_type;

      //! export the allocator type
      typedef typename BaseType:: allocator_type allocator_type;

      //! implement row_type with compressed vector
      typedef typename BaseType :: row_type row_type;

      //! increment block level counter
      enum {
        //! The number of blocklevels the matrix contains.
        blocklevel = BaseType :: blocklevel 
      };

      /** \brief Iterator for the entries of each row */
      typedef typename BaseType :: ColIterator ColIterator;

      /** \brief Iterator for the entries of each row */
      typedef typename BaseType :: ConstColIterator ConstColIterator;

      /** \brief Const iterator over the matrix rows */
      typedef typename BaseType :: RowIterator RowIterator;

      /** \brief Const iterator over the matrix rows */
      typedef typename BaseType :: ConstRowIterator ConstRowIterator;

      //! type of discrete function space 
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType SpaceType;

      //! type of block vector 
      typedef typename DiscreteFunctionType :: DofStorageType  BlockVectorType; 

      //! type of used communication manager  
      typedef CommunicationManager<SpaceType> CommunicationManagerType; 

      struct FillRowEntities
      {
        template <class GridPartImp, 
                  class EntityImp, 
                  class IndexSetImp,
                  class LocalIndicesImp>
        static void fill(const GridPartImp& gridPart,
                         const EntityImp& en,
                         const IndexSetImp& colSet,
                         LocalIndicesImp& localIndices)
        {
          const int elIndex = colSet.index( en );
          // clear set 
          localIndices.clear();
          
          // insert diagonal 
          localIndices.insert( elIndex );

          // insert neighbors 
          typedef typename GridPartImp:: IntersectionIteratorType IntersectionIteratorType; 
          IntersectionIteratorType endnit = gridPart.iend(en);
          for(IntersectionIteratorType nit = gridPart.ibegin(en); 
              nit != endnit; ++nit)
          {
            if(nit.neighbor())
            {
              const int nbIndex = colSet.index( * nit.outside() );
              localIndices.insert( nbIndex );
            }
          }
        }
      };

    private:  
      size_type nz_;

      int localRows_; 
      int localCols_;

      //! our function space, needed for communication  
      const SpaceType* space_; 

      //! communication manager 
      mutable CommunicationManagerType* comm_;

      std::vector<int> overlapRows_;

    public:
      //! constructor used by ISTLMatrixObject
      ImprovedBCRSMatrix(const SpaceType & space, 
                         CommunicationManagerType& comm,
                         size_type rows, size_type cols)
        : BaseType (rows,cols,BaseType::row_wise)
        , nz_(0)
        , space_(&space)
        , comm_(&comm)         
      {
        //std::cout << "Create Matrix with " << rows << " x " << cols << "\n";
      }

      //! constuctor used by ILU preconditioner 
      ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nz)
        : BaseType (rows,cols, BaseType::row_wise)
        , nz_(nz)
        , space_(0)
        , comm_(0)         
      {
        //std::cout << "Create Matrix with " << rows << " x " << cols << "\n";
      }
      
      //! copy constructor, needed by ISTL preconditioners 
      ImprovedBCRSMatrix(const ImprovedBCRSMatrix& org) 
        : BaseType(org) 
        , nz_(org.nz_)
        , localRows_(org.localRows_)
        , localCols_(org.localCols_)
        , space_(org.space_)
        , comm_(org.comm_)
        , overlapRows_(org.overlapRows_) 
      {}

      //! matrix multiplication for OEM solvers 
      void multOEM (const double * arg, double * dest)
      {
        RowIteratorType endi=this->end();
        for (RowIteratorType i=this->begin(); i!=endi; ++i)
        {
          int r = i.index();
          for(int k=0; k<localRows_; ++k)
          {
            int row = r * localRows_ + k; 
            dest[row] = 0.0;
          }

          ColIteratorType endj = (*i).end();
          for (ColIteratorType j=(*i).begin(); j!=endj; ++j)
          {
            for(int k=0; k<localRows_; ++k) 
            {
              int row = r * localRows_ + k; 
              for(int c=0; c<localCols_; ++c)
              {
                int col = j.index() * localCols_ + c;  
                dest[row] += (*j)[k][c] * arg[col]; 
              }
            }
          }
        }
      }

      //! setup matrix entires 
      template <class RowSpaceType, class ColSpaceType,
                class GridType,
                PartitionIteratorType pitype,
                template <class,PartitionIteratorType> class GridPartImp> 
      void setup(const RowSpaceType & rowSpace, 
                 const ColSpaceType & colSpace,
                 const GridPartImp<GridType,pitype> & gridPart,
                 bool verbose = false) 
      { 
        int size = rowSpace.indexSet().size(0);
        size = (int) size / 10;
        overlapRows_.reserve( size );
        overlapRows_.resize(0);

        assert( &rowSpace == &colSpace );
        {
          //! we need all partition iterator here  
          typedef GridPartImp<GridType,All_Partition> AllPartType; 
          typedef typename AllPartType :: template Codim<0> :: IteratorType  IteratorType;
          AllPartType allPart(const_cast<GridType&> (rowSpace.grid()));
          
          typedef typename AllPartType:: IntersectionIteratorType IntersectionIteratorType; 
          typedef typename GridType :: template Codim<0> :: Entity EntityType;
         
          typedef typename RowSpaceType :: IndexSetType RowIndexSetType;
          const RowIndexSetType& rowSet = rowSpace.indexSet();
          
          typedef typename ColSpaceType :: IndexSetType ColIndexSetType;
          const ColIndexSetType& colSet = colSpace.indexSet();

          // only works for non-hybrid grids so far 
          assert( ! rowSpace.multipleGeometryTypes() );
          
          IteratorType endit = allPart.template end<0> ();
          IteratorType it    = allPart.template begin<0> ();

          // if empty grid, do nothing
          if( it == endit ) return ;

          // initialize some values 
          localRows_ = rowSpace.baseFunctionSet(*it).numBaseFunctions();
          localCols_ = colSpace.baseFunctionSet(*it).numBaseFunctions();

          // map of indices 
          // necessary because element traversal not necessaryly is in
          // ascending order 
          std::map< int , std::set<int> > indices;
          
          for( ; it != endit; ++it)
          {
            EntityType & en = *it;
            const int elIndex = rowSet.index(en);

            // if entity is not interior, insert into overlap entries 
            if(en.partitionType() != InteriorEntity) 
              overlapRows_.push_back( elIndex );

            // set of column indices 
            std::set<int>& localIndices = indices[elIndex];

            // add all column entities to row  
            FillRowEntities::fill(allPart,en,colSet,localIndices);
          }

          typedef typename BaseType :: CreateIterator CreateIteratorType; 
          // not insert map of indices into matrix 
          CreateIteratorType endcreate = this->createend();
          for(CreateIteratorType create = this->createbegin();
              create != endcreate; ++create) 
          {
            // set of column indices 
            std::set<int>& localIndices = indices[create.index()];
            typedef typename std::set<int>::iterator iterator;
            iterator end = localIndices.end();
            // insert all indices for this row 
            for (iterator it = localIndices.begin(); it != end; ++it)
            {
              create.insert( *it );
            }
          }
        }
        std::sort(overlapRows_.begin(), overlapRows_.end());
        if(verbose)  
        {
          std::cout << "ISTLMatrix::setup: finished assembly of matrix structure! \n";
        }
      }

      //! clear Matrix, i.e. set all entires to 0
      void clear() 
      {
        {
          RowIteratorType endi=this->end();
          for (RowIteratorType i=this->begin(); i!=endi; ++i)
          {
            ColIteratorType endj = (*i).end();
            for (ColIteratorType j=(*i).begin(); j!=endj; ++j)
            {
              (*j) = 0.0;
            }
          }
        }

        // for non-interior entities set diag to 1 for ILU Preconditioner 
        const int overlap = overlapRows_.size();
        for(int i=0; i<overlap; ++i)
        {
          const int idx = overlapRows_[i]; 
          LittleBlockType& diag = this->operator[](idx)[idx];
          CompileTimeChecker<LittleBlockType :: rows == LittleBlockType :: cols > ();
          for(int k=0; k<LittleBlockType :: rows; ++k)  
          {
            diag[k][k] = 1.0;
          }
        }
      }

      //! print matrix 
      void print(std::ostream & s) const 
      {
        std::cout << "Print ISTLMatrix \n";
        RowIteratorType endi=this->end();
        for (RowIteratorType i=this->begin(); i!=endi; ++i)
        {
          ColIteratorType endj = (*i).end();
          for (ColIteratorType j=(*i).begin(); j!=endj; ++j)
          {
            s << (*j) << std::endl;
          }
        }
      }

      //! communicate block vector 
      void communicate(const BlockVectorType& arg) const 
      {
        if(comm_)
        {
          assert( space_ );
          // if serial run, just return 
          if(space_->grid().comm().size() <= 1) 
          {
            return;
          }

          // exchange data 
          DiscreteFunctionType tmp("ImprovedBCRSMatrix::communicate_tmp",*space_,arg);
          comm_->exchange( tmp );
        }
      }

      //! apply matrix: \f$ y = A(x) \f$
      void mult(const BlockVectorType& x, BlockVectorType& y) const 
      {
        // exchange data 
        communicate( x );

        // clear vector  
        y = 0;
        // multiply 
        this->umv(x,y);

        // delete non interior entries 
        deleteNonInterior(y);
      }

      //! apply scaled: \f$ y = y + \alpha A(x) \f$
      void multAdd(field_type alpha, const BlockVectorType& x, BlockVectorType& y) const 
      {
        // exchange data 
        communicate( x );

        this->usmv(alpha,x,y);

        // delete non interior entries 
        deleteNonInterior(y);
      }

    private:  
      // delete all vector entries that belong not to interior entities 
      void deleteNonInterior(BlockVectorType& y) const
      {
        // set all entries belonging to non-interior elements to zero 
        const int overlap = overlapRows_.size();
        for(int i=0; i<overlap; ++i)
        {
          y[overlapRows_[i]] = 0;
        }
      }
  };

  //! wrapper class to store perconditioner 
  //! as the interface class does not have to category 
  //! enum 
  template<class X, class Y>
  class PreconditionerWrapper : public Preconditioner<X,Y>
  {
    typedef Preconditioner<X,Y> PreconditionerInterfaceType;
    PreconditionerInterfaceType* preconder_; 
    
    //! set preconder to zero 
    PreconditionerWrapper (const PreconditionerWrapper&);
  public:
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

    enum {
      //! \brief The category the precondtioner is part of.
      category=SolverCategory::sequential};

    //! set preconder to zero 
    PreconditionerWrapper () : preconder_(0) {}
    
    //! create preconditioner of given type 
    template <class MatrixType, class PreconditionerType>
    PreconditionerWrapper(MatrixType & m, int iter, field_type relax, const PreconditionerType*) 
      : preconder_(new PreconditionerType(m,iter,relax)) {}
    
    //! create preconditioner of given type 
    template <class MatrixType, class PreconditionerType>
    PreconditionerWrapper(MatrixType & m, field_type relax, const PreconditionerType*) 
      : preconder_(new PreconditionerType(m,relax)) {}
    
    //! \copydoc Preconditioner 
    virtual void pre (X& x, Y& b) 
    {
      // apply preconditioner
      if( preconder_ ) preconder_->pre(x,b);
    }

    //! \copydoc Preconditioner 
    virtual void apply (X& v, const Y& d)
    {
      if( preconder_ ) 
      {
        // apply preconditioner
        preconder_->apply(v,d);
      }
      else 
      {
        // just copy values 
        v = d;
      }
    }

    //! \copydoc Preconditioner 
    virtual void post (X& x) {
      // apply preconditioner
      if( preconder_ ) preconder_->post(x);
    }

    // every abstract base class has a virtual destructor
    virtual ~PreconditionerWrapper () 
    {
      delete preconder_;
      preconder_ = 0;
    }
  };

  //! MatrixObject handling an istl matrix 
  template <class RowSpaceImp, class ColumnSpaceImp> 
  class ISTLMatrixObject  
  {
  public:  
    //! type of space defining row structure 
    typedef RowSpaceImp RowSpaceType;
    //! type of space defining column structure 
    typedef ColumnSpaceImp ColumnSpaceType;

  private:  
    typedef typename RowSpaceType::GridType GridType; 
    typedef typename GridType::template Codim<0>::Entity EntityType;

    enum { littleRows = RowSpaceType :: localBlockSize };
    enum { littleCols = ColumnSpaceType :: localBlockSize };
    
    typedef typename RowSpaceType :: RangeFieldType RangeFieldType;
    
    typedef FieldMatrix<RangeFieldType, littleRows, littleCols> LittleBlockType; 

    typedef BlockVectorDiscreteFunction< RowSpaceType >  DiscreteFunctionType; 
    typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType; 

  public:
    //! type of used matrix 
    typedef ImprovedBCRSMatrix< LittleBlockType , DiscreteFunctionType > MatrixType;
   
    //! type of preconditioner 
    typedef PreconditionerWrapper<BlockVectorType,BlockVectorType> PreconditionMatrixType;

    struct LocalMatrixTraits
    {
      typedef RowSpaceImp DomainSpaceType ;
      typedef ColumnSpaceImp RangeSpaceType;
      typedef double RangeFieldType;
      typedef MatrixType LocalMatrixType;
      typedef typename MatrixType:: block_type LittleBlockType;
    };

    //! LocalMatrix 
    template <class MatrixImp> 
    class LocalMatrix : public LocalMatrixDefault<LocalMatrixTraits>
    {
      typedef LocalMatrixDefault<LocalMatrixTraits> BaseType;
        
      typedef MatrixImp MatrixType;
      typedef typename MatrixType:: block_type LittleBlockType;
      
      const int rowIndex_;
      const int colIndex_;

      LittleBlockType & matrix_;
      
    public:  
      LocalMatrix(MatrixType & m,
                  const EntityType & rowEntity,
                  const RowSpaceType & rowSpace,
                  const EntityType & colEntity,
                  const ColumnSpaceType & colSpace)
        : BaseType( rowSpace, colSpace, rowEntity, colEntity )
        , rowIndex_(rowSpace.indexSet().index(rowEntity))
        , colIndex_(colSpace.indexSet().index(colEntity))
        , matrix_(m[rowIndex_][colIndex_])
      {
      }

    private: 
      LocalMatrix(const LocalMatrix&);

      // check whether given (row,col) pair is valid
      void check(int localRow, int localCol) const 
      {
        assert( localRow >= 0 );
        assert( localRow < LittleBlockType :: rows );
        assert( localCol >= 0 );
        assert( localCol < LittleBlockType :: cols );
      }
    public:
      void add(int localRow, int localCol , const double value)
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        matrix_[localRow][localCol] += value;
      }

      void set(int localRow, int localCol , const double value)
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        matrix_[localRow][localCol] = value;
      }

      double get(int localRow, int localCol ) const
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        return matrix_[localRow][localCol];
      }

      //! clear all entries belonging to local matrix 
      void clear ()
      {
        matrix_ = 0.0;
      }

      //! empty as the little matrices are already sorted
      void resort() 
      {
      }
    };

  public:
    //! type of local matrix 
    typedef LocalMatrix<MatrixType> LocalMatrixType;

  private:  
    typedef CommunicationManager<RowSpaceType> CommunicationManagerType;

    const RowSpaceType & rowSpace_;
    const ColumnSpaceType & colSpace_;

    int size_;

    int sequence_;

    mutable MatrixType* matrix_;
    mutable PreconditionMatrixType* preconder_;

    CommunicationManagerType comm_;

    int numIterations_; 
    double relaxFactor_; 
      
    enum PreConder_Id { none  = 0 , // no preconditioner 
                        ssor  = 1 , // SSOR preconditioner 
                        sor   = 2 , // SOR preconditioner 
                        ilu_0 = 3 , // ILU-0 preconditioner 
                        ilu_n = 4 , // ILU-n preconditioner 
                        gauss_seidel= 5 , // Gauss-Seidel preconditioner 
                        jacobi = 6  // Jacobi preconditioner 
    };
    
    PreConder_Id preconditioning_;

    // prohibit copy constructor 
    ISTLMatrixObject(const ISTLMatrixObject&); 
  public:  
    //! constructor 
    //! \param rowSpace space defining row structure 
    //! \param colSpace space defining column structure 
    //! \param paramfile parameter file to read variables 
    //!         - Preconditioning: {0,1,2,3,4,5,6} put -1 to get info
    //!         - Pre-iteration: number of iteration of preconditioner
    //!         - Pre-relaxation: relaxation factor   
    ISTLMatrixObject(const RowSpaceType & rowSpace,
                     const ColumnSpaceType & colSpace,
                     const std::string& paramfile)
      : rowSpace_(rowSpace)
      , colSpace_(colSpace)
      , size_(-1)
      , sequence_(-1)
      , matrix_(0)
      , preconder_(0)
      , comm_(rowSpace_)
      , numIterations_(5)
      , relaxFactor_(1.1)
      , preconditioning_(none)
    {
      if(paramfile != "")
      {
        int preCon = 0;
        readParameter(paramfile,"Preconditioning",preCon);
        if( preCon >= 0 && preCon <= 6) 
          preconditioning_ = (PreConder_Id) preCon;
        else 
          preConErrorMsg(preCon);
        
        readParameter(paramfile,"Pre-iteration",numIterations_);
        readParameter(paramfile,"Pre-relaxation",relaxFactor_);
      }

      // only ILU-0 works savely in parallel
      if(rowSpace_.grid().comm().size() > 1)
      {
        if( (preconditioning_ != none) && (preconditioning_ != ilu_0) )
        {
          std::cerr << "ERROR: Only Preconditioner ILU-0 works in parallel! " << std::endl;
          abort();
        }
      }
      
      assert( rowSpace_.indexSet().size(0) ==
              colSpace_.indexSet().size(0) );
    }

    //! destructor 
    ~ISTLMatrixObject() 
    {
      delete preconder_;
      delete matrix_;
    }

    //! return reference to system matrix 
    MatrixType & matrix() const 
    { 
      assert( matrix_ );
      return *matrix_; 
    }
    
    //! return true, because in case of no preconditioning we have empty
    //! preconditioner 
    bool hasPcMatrix () const { return true; }

    //! return reference to preconditioner
    const PreconditionMatrixType& pcMatrix () const  
    { 
      if( !preconder_ )
      {
        preconder_ = createPreconditioner();
      }
      return *preconder_; 
    }

    //! set all matrix entries to zero 
    void clear()
    {
      matrix().clear();
    }

    //! reserve matrix with right size 
    void reserve(bool verbose = false) 
    {
      // if grid sequence number changed, rebuild matrix 
      if(sequence_ != rowSpace_.sequence())
      {
        delete matrix_; matrix_ = 0;
        delete preconder_; preconder_ = 0;
        size_ = rowSpace_.indexSet().size(0);

        matrix_ = new MatrixType(rowSpace_, comm_, size_, colSpace_.indexSet().size(0));
        matrix().setup(rowSpace_,colSpace_,rowSpace_.gridPart(),verbose);

        sequence_ = rowSpace_.sequence();
      }
    }

    //! mult method of matrix object used by oem solver
    void multOEM(const double * arg, double * dest) const
    {
      matrix().multOEM(arg,dest);
    }

    //! resort row numbering in matrix to have ascending numbering 
    void resort()
    {
    }

    //! create precondition matrix does nothing because preconditioner is
    //! created only when requested 
    void createPreconditionMatrix() 
    {
    }

    //! print matrix 
    void print(std::ostream & s) const 
    { 
      matrix().print(std::cout);
    }

  private:  
    PreconditionMatrixType* createPreconditioner() const
    {
      // no preconditioner 
      if( preconditioning_ == none )
      {
        return new PreconditionMatrixType();
      }
      // SSOR 
      else if( preconditioning_ == ssor )
      {
        typedef SeqSSOR<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      // SOR 
      else if(preconditioning_ == sor )
      {
        typedef SeqSOR<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      // ILU-0 
      else if(preconditioning_ == ilu_0)
      {
        typedef SeqILU0<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), relaxFactor_, (PreconditionerType*)0);
      }
      // ILU-n
      else if(preconditioning_ == ilu_n)
      {
        typedef SeqILUn<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      // Gauss-Seidel
      else if(preconditioning_ == gauss_seidel)
      {
        typedef SeqGS<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      // Jacobi 
      else if(preconditioning_ == jacobi)
      {
        typedef SeqJac<MatrixType,BlockVectorType,BlockVectorType> PreconditionerType;
        return new PreconditionMatrixType(matrix(), numIterations_ , relaxFactor_, (PreconditionerType*)0);
      }
      else 
      {
        preConErrorMsg(preconditioning_);
      }
      return 0;
    }

    void preConErrorMsg(int preCon) const 
    {
      std::cerr << "ERROR: Wrong precoditioning number (p = " << preCon;
      std::cerr <<") in ISTLMatrixObject! \n";
      std::cerr <<"Valid values are: \n";
      std::cerr <<"0 == no \n";
      std::cerr <<"1 == SSOR \n";
      std::cerr <<"2 == SOR \n";
      std::cerr <<"3 == ILU-0 \n";
      std::cerr <<"4 == ILU-n \n";
      std::cerr <<"5 == Gauss-Seidel \n";
      std::cerr <<"6 == Jacobi \n";
      assert(false);
      exit(1);
    }
  };

} // end namespace Dune 
#endif

#endif
