#ifndef DUNE_ISTLMATRIXWRAPPER_HH
#define DUNE_ISTLMATRIXWRAPPER_HH

#if HAVE_DUNE_ISTL

//- system includes 
#include <vector> 

//- Dune common includes 
#include <dune/common/exceptions.hh>

//- Dune istl includes 
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>

//- Dune fem includes 
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/matrix/preconditionerwrapper.hh>
#include <dune/fem/io/parameter.hh>

namespace Dune
{ 
  namespace Fem
  {

  ///////////////////////////////////////////////////////
  // --ISTLMatrixHandle
  //////////////////////////////////////////////////////
  template <class LittleBlockType, 
            class RowDiscreteFunctionImp, 
            class ColDiscreteFunctionImp = RowDiscreteFunctionImp> 
  class ImprovedBCRSMatrix : public BCRSMatrix<LittleBlockType> 
  {
    public:
      typedef RowDiscreteFunctionImp RowDiscreteFunctionType;
      typedef ColDiscreteFunctionImp ColDiscreteFunctionType;
      
      typedef BCRSMatrix<LittleBlockType> BaseType; 
      typedef typename BaseType :: RowIterator RowIteratorType ;
      typedef typename BaseType :: ColIterator ColIteratorType ;

      typedef ImprovedBCRSMatrix< LittleBlockType,
              RowDiscreteFunctionImp, ColDiscreteFunctionImp > ThisType;

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
      typedef typename ColDiscreteFunctionType :: DiscreteFunctionSpaceType RangeSpaceType;

      //! type of row block vector 
      typedef typename RowDiscreteFunctionType :: DofStorageType  RowBlockVectorType; 

      //! type of column block vector 
      typedef typename ColDiscreteFunctionType :: DofStorageType  ColBlockVectorType; 

      //! type of communication object 
      typedef typename RangeSpaceType :: GridType :: Traits :: CollectiveCommunication   CollectiveCommunictionType ;

      typedef typename BaseType :: BuildMode BuildMode ;

    public:
      //! constructor used by ISTLMatrixObject
      ImprovedBCRSMatrix(size_type rows, size_type cols) 
        : BaseType (rows,cols, BaseType :: row_wise)
      {
      }

      //! constuctor used by ILU preconditioner 
      ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nz) 
        : BaseType (rows,cols, BaseType :: row_wise)
      {
      }
      
      //! copy constructor, needed by ISTL preconditioners 
      ImprovedBCRSMatrix( ) 
        : BaseType ()
      {}

      //! copy constructor, needed by ISTL preconditioners 
      ImprovedBCRSMatrix(const ImprovedBCRSMatrix& org) 
        : BaseType(org) 
      {}

      //! setup matrix entires 
      template <class RowMapperType, class ColMapperType,
                class StencilCreatorImp> 
      void setup(const RangeSpaceType& colSpace, 
                 const RowMapperType & rowMapper, 
                 const ColMapperType & colMapper,
                 const StencilCreatorImp& stencil, 
                 bool verbose = false) 
      { 
        // if empty grid, do nothing
        if( colSpace.begin() == colSpace.end() ) return ;
        
        {
          // map of indices 
          // necessary because element traversal not necessaryly is in
          // ascending order 
          std::map< int , std::set<int> > indices;

          // build matrix entries
          stencil.setup(colSpace, rowMapper, colMapper, indices , (ColDiscreteFunctionType*) 0);

          // insert entries 
          createEntries( indices );
        }

        // in verbose mode some output 
        if(verbose)  
        {
          std::cout << "ISTLMatrix::setup: finished assembly of matrix structure! \n";
        }
      }

      void createEntries(std::map<int , std::set<int> >& indices) 
      {
        // type of create interator 
        typedef typename BaseType :: CreateIterator CreateIteratorType; 
        // not insert map of indices into matrix 
        CreateIteratorType endcreate = this->createend();
        for(CreateIteratorType create = this->createbegin();
            create != endcreate; ++create) 
        {
          // set of column indices 
          std::set<int>& localIndices = indices[ create.index() ];
          typedef typename std::set<int>::iterator iterator;
          iterator end = localIndices.end();
          // insert all indices for this row 
          for (iterator it = localIndices.begin(); it != end; ++it)
          {
            create.insert( *it );
          }
        }
      }

      //! clear Matrix, i.e. set all entires to 0
      void clear() 
      {
        RowIteratorType endi=this->end();
        for (RowIteratorType i=this->begin(); i!=endi; ++i)
        {
          ColIteratorType endj = (*i).end();
          for (ColIteratorType j=(*i).begin(); j!=endj; ++j)
          {
            (*j) = 0;
          }
        }
      }

      //! setup like the old matrix but remove rows with hanging nodes 
      template <class HangingNodesType> 
      void setup(ThisType& oldMatrix,
                 const HangingNodesType& hangingNodes) 
      {
        // necessary because element traversal not necessaryly is in
        // ascending order 
        typedef std::set< std::pair<int, block_type> > LocalEntryType;
        typedef std::map< int , LocalEntryType > EntriesType;
        EntriesType entries;

        {
          // map of indices 
          std::map< int , std::set<int> > indices;
          // not insert map of indices into matrix 
          RowIteratorType rowend  = oldMatrix.end();
          for(RowIteratorType it  = oldMatrix.begin(); it != rowend; ++it)
          {
            const int row = it.index();
            std::set< int >& localIndices = indices[ row ];

            if( hangingNodes.isHangingNode( row ) )
            {
              // insert columns into other columns 
              typedef typename HangingNodesType :: ColumnVectorType ColumnVectorType;
              const ColumnVectorType& cols = hangingNodes.associatedDofs( row );
              const size_t colSize = cols.size();
              for(size_t i=0; i<colSize; ++i) 
              {
                assert( ! hangingNodes.isHangingNode( cols[i].first ) );

                // get local indices of col
                std::set< int >& localColIndices = indices[ cols[i].first ];
                LocalEntryType& localEntry = entries[  cols[i].first ];

                // copy from old matrix 
                ColIteratorType endj = (*it).end();
                for (ColIteratorType j= (*it).begin(); j!=endj; ++j)
                {
                  localColIndices.insert( j.index () );
                  localEntry.insert( std::make_pair( j.index(), (cols[i].second * (*j)) ));
                }
              }

              // insert diagonal and hanging columns 
              localIndices.insert( row );
              for(size_t i=0; i<colSize; ++i) 
              {
                localIndices.insert( cols[i].first );
              }
            }
            else 
            {
              // copy from old matrix 
              ColIteratorType endj = (*it).end();
              for (ColIteratorType j= (*it).begin(); j!=endj; ++j)
              {
                localIndices.insert( j.index () );
              }
            }
          }

          // create matrix from entry map 
          createEntries( indices );

        } // end create, matrix is on delete of create iterator 

        {
          // not insert map of indices into matrix 
          RowIteratorType rowit  = oldMatrix.begin();

          RowIteratorType endcreate = this->end();
          for(RowIteratorType create = this->begin();
              create != endcreate; ++create, ++rowit ) 
          {
            assert( rowit != oldMatrix.end() );

            const int row = create.index();
            if( hangingNodes.isHangingNode( row ) )
            {
              typedef typename HangingNodesType :: ColumnVectorType ColumnVectorType;
              const ColumnVectorType& cols = hangingNodes.associatedDofs( row );

              std::map< const int , block_type > colMap; 
              // only working for block size 1 ath the moment 
              assert( block_type :: rows == 1 );
              // insert columns into map 
              const size_t colSize = cols.size();
              for( size_t i=0; i<colSize; ++i) 
              {
                colMap[ cols[i].first ] = -cols[i].second;
              }
              // insert diagonal into map
              colMap[ row ] = 1;

              ColIteratorType endj = (*create).end();
              for (ColIteratorType j= (*create).begin(); j!=endj; ++j)
              {
                assert( colMap.find( j.index() ) != colMap.end() );
                (*j) = colMap[ j.index() ];
              }
            }
            // if entries are equal, just copy 
            else if ( entries.find( row ) == entries.end() ) 
            {
              ColIteratorType colit = (*rowit).begin();
              ColIteratorType endj = (*create).end();
              for (ColIteratorType j= (*create).begin(); j!=endj; ++j, ++colit )
              {
                assert( colit != (*rowit).end() );
                (*j) = (*colit);
              }
            }
            else 
            {
              typedef std::map< int , block_type > ColMapType;
              ColMapType oldCols;

              {
                ColIteratorType colend = (*rowit).end(); 
                for(ColIteratorType colit = (*rowit).begin(); colit !=
                    colend; ++colit)
                {
                  oldCols[ colit.index() ] = 0; 
                }
              }

              typedef typename EntriesType :: iterator Entryiterator ;
              Entryiterator entry = entries.find( row );
              assert( entry  != entries.end ());

              {
                typedef typename LocalEntryType :: iterator iterator;
                iterator endcol = (*entry).second.end();
                for( iterator co = (*entry).second.begin(); co != endcol; ++co) 
                {
                  oldCols[ (*co).first ] = 0;
                }
              }

              {
                ColIteratorType colend = (*rowit).end(); 
                for(ColIteratorType colit = (*rowit).begin(); colit !=
                    colend; ++colit)
                {
                  oldCols[ colit.index() ] += (*colit); 
                }
              }

              {
                typedef typename LocalEntryType :: iterator iterator;
                iterator endcol = (*entry).second.end();
                for( iterator co = (*entry).second.begin(); co != endcol; ++co) 
                {
                  oldCols[ (*co).first ] += (*co).second;   
                }
              }

              ColIteratorType endj = (*create).end();
              for (ColIteratorType j= (*create).begin(); j!=endj; ++j )
              {
                typedef typename ColMapType :: iterator iterator;
                iterator colEntry = oldCols.find( j.index() );
                if( colEntry != oldCols.end() ) 
                {
                  (*j) = (*colEntry).second;
                }
                else 
                {
                  abort();
                }
              }
            }
          }
        } // end create 
      }

      //! extract diagonal of matrix to block vector 
      void extractDiagonal( ColBlockVectorType& diag ) const
      { 
        ConstRowIterator endi = this->end();
        for (ConstRowIterator i = this->begin(); i!=endi; ++i)
        {
          const size_t row = i.index();
          const ConstColIterator endj = (*i).end();
          for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          {
            if( j.index() == row ) 
            {
              const LittleBlockType& block = (*j);
              enum { blockSize = LittleBlockType :: rows };
              for( int l=0; l<blockSize; ++l )
              {
                diag[ row ][ l ] = block[ l ][ l ]; 
              }
            }
          }
        }
      }

      //! print matrix 
      void print(std::ostream & s) const 
      {
        std::cout << "Print ISTLMatrix \n";
        ConstRowIterator endi=this->end();
        for (ConstRowIterator i=this->begin(); i!=endi; ++i)
        {
          ConstColIterator endj = (*i).end();
          for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
          {
            s << (*j) << std::endl;
          }
        }
      }
  };

  template <class RowSpaceImp, class ColSpaceImp, class TraitsImp> 
  class ISTLMatrixObject;
  
  template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
  struct ISTLMatrixTraits
  {
    typedef RowSpaceImp DomainSpaceType; 
    typedef ColSpaceImp RangeSpaceType; 
    typedef ISTLMatrixTraits<DomainSpaceType,RangeSpaceType> ThisType;  
    
    template <class OperatorTraits>
    struct MatrixObject 
    {
      typedef ISTLMatrixObject<DomainSpaceType,RangeSpaceType,OperatorTraits> MatrixObjectType; 
    };
  };

  //! MatrixObject handling an istl matrix 
  template <class DomainSpaceImp, class RangeSpaceImp, class TraitsImp> 
  class ISTLMatrixObject
  {
  public:  
    //! type of traits 
    typedef TraitsImp Traits;

    //! type of stencil defined by operator
    typedef typename Traits :: StencilType StencilType;
    
    //! type of space defining row structure 
    typedef DomainSpaceImp DomainSpaceType; 
    //typedef typename Traits :: DomainSpaceType DomainSpaceType;
    //! type of space defining column structure 
    typedef RangeSpaceImp RangeSpaceType; 
    //typedef typename Traits :: RangeSpaceType RangeSpaceType;

    //! type of this pointer 
    typedef ISTLMatrixObject<DomainSpaceType,RangeSpaceType,Traits> ThisType;


  protected:  
    typedef typename DomainSpaceType::GridType GridType; 

    typedef typename RangeSpaceType :: EntityType ColumnEntityType ;
    typedef typename DomainSpaceType :: EntityType RowEntityType ;

    enum { littleRows = DomainSpaceType :: localBlockSize };
    enum { littleCols = RangeSpaceType :: localBlockSize };
    
    typedef typename DomainSpaceType :: RangeFieldType RangeFieldType;
    
    typedef FieldMatrix<RangeFieldType, littleRows, littleCols> LittleBlockType; 

  public:
    typedef ISTLBlockVectorDiscreteFunction< DomainSpaceType >     RowDiscreteFunctionType; 
    typedef typename RowDiscreteFunctionType :: LeakPointerType  RowLeakPointerType;
    typedef ISTLBlockVectorDiscreteFunction< RangeSpaceType >  ColumnDiscreteFunctionType; 
    typedef typename ColumnDiscreteFunctionType :: LeakPointerType  ColumnLeakPointerType;
    
  protected:  
    typedef typename RowDiscreteFunctionType :: DofStorageType    RowBlockVectorType; 
    typedef typename ColumnDiscreteFunctionType :: DofStorageType ColumnBlockVectorType; 

    typedef typename DomainSpaceType :: BlockMapperType RowMapperType; 
    typedef typename RangeSpaceType :: BlockMapperType ColMapperType; 

  public:
    //! type of used matrix 
    typedef ImprovedBCRSMatrix< LittleBlockType , 
                                RowDiscreteFunctionType , 
                                ColumnDiscreteFunctionType > MatrixType;
    typedef typename Traits :: template Adapter < MatrixType > ::  MatrixAdapterType MatrixAdapterType;
    // get preconditioner type from MatrixAdapterType
    typedef ThisType PreconditionMatrixType;
    typedef typename MatrixAdapterType :: ParallelScalarProductType ParallelScalarProductType;
   
    template <class MatrixObjectImp> 
    class LocalMatrix; 

    struct LocalMatrixTraits
    {
      typedef DomainSpaceImp DomainSpaceType ;
      typedef RangeSpaceImp  RangeSpaceType;
      typedef typename DomainSpaceType :: RangeFieldType RangeFieldType;
      typedef LocalMatrix<ThisType> LocalMatrixType;
      typedef typename MatrixType:: block_type LittleBlockType;
    };

    //! LocalMatrix 
    template <class MatrixObjectImp> 
    class LocalMatrix : public LocalMatrixDefault<LocalMatrixTraits>
    {
    public:  
      //! type of base class 
      typedef LocalMatrixDefault<LocalMatrixTraits> BaseType;

      //! type of matrix object 
      typedef MatrixObjectImp MatrixObjectType;
      //! type of matrix 
      typedef typename MatrixObjectImp :: MatrixType MatrixType;
      //! type of little blocks 
      typedef typename MatrixType:: block_type LittleBlockType;
      //! type of entries of little blocks 
      typedef typename DomainSpaceType :: RangeFieldType DofType;

      typedef typename MatrixType::row_type RowType; 

      //! type of row mapper 
      typedef typename MatrixObjectType :: RowMapperType RowMapperType;
      //! type of col mapper 
      typedef typename MatrixObjectType :: ColMapperType ColMapperType;
        
    private:
      // special mapper omiting block size 
      const RowMapperType& rowMapper_;
      const ColMapperType& colMapper_;
      
      // number of local matrices 
      int numRows_;
      int numCols_;

      // vector with pointers to local matrices 
      typedef std::vector<LittleBlockType *> LittleMatrixRowStorageType ;
      std::vector< LittleMatrixRowStorageType > matrices_;

      // matrix to build 
      const MatrixObjectType& matrixObj_;

      // type of actual geometry 
      GeometryType geomType_;
      
    public:  
      LocalMatrix(const MatrixObjectType & mObj,
                  const DomainSpaceType & rowSpace,
                  const RangeSpaceType & colSpace)
        : BaseType( rowSpace, colSpace )
        , rowMapper_(mObj.rowMapper())
        , colMapper_(mObj.colMapper())
        , numRows_( rowMapper_.maxNumDofs() )
        , numCols_( colMapper_.maxNumDofs() )
        , matrixObj_(mObj)
        , geomType_(GeometryType::simplex,0)
      {
      }

      void init(const RowEntityType & rowEntity,
                const ColumnEntityType & colEntity)
      {
        if( geomType_ != rowEntity.type() ) 
        {
          // initialize base functions sets 
          BaseType :: init ( rowEntity , colEntity );

          geomType_ = rowEntity.type();
          numRows_  = rowMapper_.numDofs(rowEntity);
          numCols_  = colMapper_.numDofs(colEntity);
          matrices_.resize( numRows_ );
        }

        MatrixType& matrix = matrixObj_.matrix();
        typedef typename RowMapperType :: DofMapIteratorType RowMapIteratorType ;
        typedef typename ColMapperType :: DofMapIteratorType ColMapIteratorType ;

        const RowMapIteratorType endrow = rowMapper_.end( rowEntity );
        for( RowMapIteratorType row = rowMapper_.begin( rowEntity );
             row != endrow; ++row ) 
        {
          LittleMatrixRowStorageType& localMatRow = matrices_[ row.local() ];
          localMatRow.resize( numCols_ );

          // get row 
          RowType& matRow = matrix[ row.global() ];

          const ColMapIteratorType endcol = colMapper_.end( colEntity );
          for( ColMapIteratorType col = colMapper_.begin( colEntity );
               col != endcol; ++col ) 
          {
            assert( matrix.exists( row.global(), col.global() ) );
            localMatRow[ col.local() ] = &matRow[ col.global() ];
          }
        }
      }

      LocalMatrix(const LocalMatrix& org) 
        : BaseType( org )
        , rowMapper_(org.rowMapper_)
        , colMapper_(org.colMapper_)
        , numRows_( org.numRows_ )
        , numCols_( org.numCols_ ) 
        , matrices_(org.matrices_)
        , matrixObj_(org.matrixObj_)
        , geomType_(org.geomType_)
      {
      }

    private: 
      // check whether given (row,col) pair is valid
      void check(int localRow, int localCol) const 
      {
        const size_t row = (int) localRow / littleRows;
        const size_t col = (int) localCol / littleCols;
        const int lRow = localRow%littleRows;
        const int lCol = localCol%littleCols;
        assert( row < matrices_.size() ) ;
        assert( col < matrices_[row].size() ); 
        assert( lRow < littleRows );
        assert( lCol < littleCols ); 
      }

      DofType& getValue(const int localRow, const int localCol) 
      {
        const int row = (int) localRow / littleRows;
        const int col = (int) localCol / littleCols;
        const int lRow = localRow%littleRows;
        const int lCol = localCol%littleCols;
        return (*matrices_[row][col])[lRow][lCol];
      }

    public:
      const DofType get(const int localRow, const int localCol) const 
      {
        const int row = (int) localRow / littleRows;
        const int col = (int) localCol / littleCols;
        const int lRow = localRow%littleRows;
        const int lCol = localCol%littleCols;
        return (*matrices_[row][col])[lRow][lCol];
      }

      void scale (const DofType& scalar) 
      {
        for(size_t i=0; i<matrices_.size(); ++i)
          for(size_t j=0; j<matrices_[i].size(); ++j)
            (*matrices_[i][j]) *= scalar;
      }

      void add(const int localRow, const int localCol , const DofType value)
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        getValue(localRow,localCol) += value;
      }

      void set(const int localRow, const int localCol , const DofType value)
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        getValue(localRow,localCol) = value;
      }

      //! make unit row (all zero, diagonal entry 1.0 )
      void unitRow(const int localRow) 
      {
        const int row = (int) localRow / littleRows;
        const int lRow = localRow%littleRows;

        // clear row 
        doClearRow( row, lRow );

        // set diagonal entry to 1 
        (*matrices_[row][row])[lRow][lRow] = 1;
      }

      //! clear all entries belonging to local matrix 
      void clear ()
      {
        for(int i=0; i<matrices_.size(); ++i)
          for(int j=0; j<matrices_[i].size(); ++j)
            (*matrices_[i][j]) = (DofType) 0;
      }

      //! set matrix row to zero
      void clearRow ( const int localRow )
      {
        const int row = (int) localRow / littleRows;
        const int lRow = localRow%littleRows;

        // clear the row 
        doClearRow( row, lRow );
      }

      //! empty as the little matrices are already sorted
      void resort ()
      {}

  protected:  
      //! set matrix row to zero
      void doClearRow ( const int row, const int lRow )
      {
        // get number of columns  
        const int col = this->columns();
        for(int localCol=0; localCol<col; ++localCol) 
        {
          const int col = (int) localCol / littleCols;
          const int lCol = localCol%littleCols;
          (*matrices_[row][col])[lRow][lCol] = 0;
        }
      }

    }; // end of class LocalMatrix 

  public:
    //! type of local matrix 
    typedef LocalMatrix<ThisType> ObjectType;
    typedef ThisType LocalMatrixFactoryType;
    typedef Fem :: ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
    //! type of local matrix 
    typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

  protected:  
    const DomainSpaceType & domainSpace_;
    const RangeSpaceType & rangeSpace_;

    // sepcial row mapper 
    RowMapperType& rowMapper_;
    // special col mapper 
    ColMapperType& colMapper_;

    int size_;

    int sequence_;

    mutable MatrixType* matrix_;

    ParallelScalarProductType scp_;

    int numIterations_; 
    double relaxFactor_; 
      
    enum ISTLPreConder_Id { none  = 0 ,      // no preconditioner 
                            ssor  = 1 ,      // SSOR preconditioner 
                            sor   = 2 ,      // SOR preconditioner 
                            ilu_0 = 3 ,      // ILU-0 preconditioner 
                            ilu_n = 4 ,      // ILU-n preconditioner 
                            gauss_seidel= 5, // Gauss-Seidel preconditioner 
                            jacobi = 6,      // Jacobi preconditioner 
                            amg_ilu_0 = 7,   // AMG with ILU-0 smoother 
                            amg_ilu_n = 8,   // AMG with ILU-n smoother  
                            amg_jacobi = 9   // AMG with Jacobi smoother  
    };
    
    ISTLPreConder_Id preconditioning_;

    mutable LocalMatrixStackType localMatrixStack_;

    mutable MatrixAdapterType* matrixAdap_;
    mutable RowBlockVectorType* Arg_;
    mutable ColumnBlockVectorType* Dest_;

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
    ISTLMatrixObject ( const DomainSpaceType &rowSpace,
                       const RangeSpaceType &colSpace,
                       const std :: string &paramfile = "" )
      : domainSpace_(rowSpace)
      , rangeSpace_(colSpace)
      // create scp to have at least one instance 
      // otherwise instance will be deleted during setup
      // get new mappers with number of dofs without considerung block size 
      , rowMapper_( rowSpace.blockMapper() )
      , colMapper_( colSpace.blockMapper() )
      , size_(-1)
      , sequence_(-1)
      , matrix_(0)
      , scp_(rangeSpace())
      , numIterations_(5)
      , relaxFactor_(1.1)
      , preconditioning_(none)
      , localMatrixStack_( *this )
      , matrixAdap_(0)
      , Arg_(0)
      , Dest_(0)
    {
      int preCon = 0;
      if(paramfile != "")
      {
        DUNE_THROW(InvalidStateException,"ISTLMatrixObject: old parameter method disabled");
      }
      else 
      {
        static const std::string preConTable[]
          = { "none", "ssor", "sor", "ilu-0", "ilu-n", "gauss-seidel", "jacobi",
            "amg-ilu-0", "amg-ilu-n", "amg-jacobi" };
        preCon         = Parameter::getEnum( "istl.preconditioning.method", preConTable, preCon );
        numIterations_ = Parameter::getValue( "istl.preconditioning.iterations", numIterations_ );
        relaxFactor_   = Parameter::getValue( "istl.preconditioning.relaxation", relaxFactor_ );
      }

      if( preCon >= 0 && preCon <= 9) 
        preconditioning_ = (ISTLPreConder_Id) preCon;
      else 
        preConErrorMsg(preCon);

      assert( rowMapper_.size() == colMapper_.size() );
    }

    /** \copydoc Fem::Operator::assembled */
    static const bool assembled = true ;

    const ThisType& systemMatrix() const { return *this; }
  public:  
    //! destructor 
    ~ISTLMatrixObject() 
    {
      removeObj( true );
    }

    //! return reference to system matrix 
    MatrixType & matrix() const 
    { 
      assert( matrix_ );
      return *matrix_; 
    }

    void printTexInfo(std::ostream& out) const
    {
      out << "ISTL MatrixObj: ";
      out << " preconditioner = " << preconditionName() ;
      out  << "\\\\ \n";
    }

    //! return matrix adapter object  
    std::string preconditionName() const 
    { 
      std::stringstream tmp ; 
      // no preconditioner 
      switch (preconditioning_) 
      {
        case ssor : tmp << "SSOR"; break;
        case sor  : tmp << "SOR"; break;
        case ilu_0: tmp << "ILU-0"; break;
        case ilu_n: tmp << "ILU-n"; break;
        case gauss_seidel : tmp << "Gauss-Seidel"; break;
        case jacobi: tmp << "Jacobi"; break; 
        default: tmp << "None"; break;
      }

      if( preconditioning_ != ilu_0 ) 
      {
        tmp << " n=" << numIterations_;
      }
      tmp << " relax=" << relaxFactor_ ;
      return tmp.str();
    }
    
    template <class PreconditionerType> 
    MatrixAdapterType 
    createMatrixAdapter(const PreconditionerType* preconditioning,
                        size_t numIterations) const 
    {
      typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
      PreConType preconAdapter(matrix(), numIterations, relaxFactor_, preconditioning );
      return MatrixAdapterType(matrix(), domainSpace(), rangeSpace(), preconAdapter );
    }

    template <class PreconditionerType> 
    MatrixAdapterType 
    createAMGMatrixAdapter(const PreconditionerType* preconditioning,
                           size_t numIterations) const 
    {
      typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
      PreConType preconAdapter(matrix(), numIterations, relaxFactor_, preconditioning, domainSpace().grid().comm() );
      return MatrixAdapterType(matrix(), domainSpace(), rangeSpace(), preconAdapter );
    }

    //! return matrix adapter object  
    const MatrixAdapterType& matrixAdapter() const 
    {
      if( matrixAdap_ == 0 ) 
        matrixAdap_ = new MatrixAdapterType( matrixAdapterObject() );
      return *matrixAdap_;
    }

  protected:  
    MatrixAdapterType matrixAdapterObject() const 
    { 
#ifndef DISABLE_ISTL_PRECONDITIONING
      const size_t procs = domainSpace().grid().comm().size();

      typedef typename MatrixType :: BaseType ISTLMatrixType ;
      typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
      // no preconditioner 
      if( preconditioning_ == none )
      {
        return MatrixAdapterType(matrix(), domainSpace(),rangeSpace(), PreConType() );
      }
      // SSOR 
      else if( preconditioning_ == ssor )
      {
        if( procs > 1 ) 
          DUNE_THROW(InvalidStateException,"ISTL::SeqSSOR not working in parallel computations");

        typedef SeqSSOR<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return createMatrixAdapter( (PreconditionerType*)0, numIterations_ );
      }
      // SOR 
      else if(preconditioning_ == sor )
      {
        if( procs > 1 ) 
          DUNE_THROW(InvalidStateException,"ISTL::SeqSOR not working in parallel computations");

        typedef SeqSOR<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return createMatrixAdapter( (PreconditionerType*)0, numIterations_ );
      }
      // ILU-0 
      else if(preconditioning_ == ilu_0)
      {
        if( procs > 1 ) 
          DUNE_THROW(InvalidStateException,"ISTL::SeqILU0 not working in parallel computations");

        typedef FemSeqILU0<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return createMatrixAdapter( (PreconditionerType*)0, numIterations_ );
      }
      // ILU-n
      else if(preconditioning_ == ilu_n)
      {
        if( procs > 1 ) 
          DUNE_THROW(InvalidStateException,"ISTL::SeqILUn not working in parallel computations");

        typedef SeqILUn<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return createMatrixAdapter( (PreconditionerType*)0, numIterations_ );
      }
      // Gauss-Seidel
      else if(preconditioning_ == gauss_seidel)
      {
        if( procs > 1 ) 
          DUNE_THROW(InvalidStateException,"ISTL::SeqGS not working in parallel computations");

        typedef SeqGS<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return createMatrixAdapter( (PreconditionerType*)0, numIterations_ );
      }
      // Jacobi 
      else if(preconditioning_ == jacobi)
      {
        if( numIterations_ == 1 ) // diagonal preconditioning 
        {
          typedef FemDiagonalPreconditioner< ThisType, RowBlockVectorType, ColumnBlockVectorType > PreconditionerType;
          typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
          PreConType preconAdapter( matrix(), new PreconditionerType( *this ) );
          return MatrixAdapterType( matrix(), domainSpace(), rangeSpace(), preconAdapter );
        }
        else if ( procs == 1 )
        {
          typedef SeqJac<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (PreconditionerType*)0, numIterations_ );
        }
        else 
        {
          DUNE_THROW(InvalidStateException,"ISTL::SeqJac(Jacobi) only working with istl.preconditioning.iterations: 1 in parallel computations");
        }
      }
      // AMG ILU-0  
      else if(preconditioning_ == amg_ilu_0)
      {
        // use original SeqILU0 because of some AMG traits classes.
        typedef SeqILU0<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return createAMGMatrixAdapter( (PreconditionerType*)0, numIterations_ );
      }
      // AMG ILU-n  
      else if(preconditioning_ == amg_ilu_n)
      {
        typedef SeqILUn<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return createAMGMatrixAdapter( (PreconditionerType*)0, numIterations_ );
      }
      // AMG Jacobi   
      else if(preconditioning_ == amg_jacobi)
      {
        typedef SeqJac<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
        return createAMGMatrixAdapter( (PreconditionerType*)0, numIterations_ );
      }
      else 
      {
        preConErrorMsg(preconditioning_);
      }
#endif

      return MatrixAdapterType(matrix(), domainSpace(), rangeSpace(), PreConType() );
    }
    
  public:  
    //! return true, because in case of no preconditioning we have empty
    //! preconditioner (used by OEM methods)
    bool hasPreconditionMatrix() const { return (preconditioning_ != none); }

    //! return reference to preconditioner object (used by OEM methods)
    const PreconditionMatrixType& preconditionMatrix() const { return *this; }

    //! set all matrix entries to zero 
    void clear()
    {
      matrix().clear();
      // clean matrix adapter and other helper classes 
      removeObj( false );
    }

    //! reserve matrix with right size 
    void reserve(bool verbose = false) 
    {
      // if grid sequence number changed, rebuild matrix 
      if(sequence_ != domainSpace().sequence())
      {
        removeObj( true );

        StencilType stencil; 
        matrix_ = new MatrixType(rowMapper_.size(),colMapper_.size());
        matrix().setup(rangeSpace(),rowMapper(),colMapper(),stencil,verbose);

        sequence_ = domainSpace().sequence();
      }
    }

    //! setup new matrix with hanging nodes taken into account 
    template <class HangingNodesType> 
    void changeHangingNodes(const HangingNodesType& hangingNodes) 
    {
      // create new matrix 
      MatrixType* newMatrix = new MatrixType(rowMapper_.size(), colMapper_.size());

      // setup with hanging rows 
      newMatrix->setup( *matrix_ , hangingNodes );

      // remove old matrix 
      removeObj( true );
      // store new matrix 
      matrix_ = newMatrix;
    }

    //! extract diagonal entries of the matrix to a discrete function of type
    //! BlockVectorDiscreteFunction 
    void extractDiagonal( ColumnDiscreteFunctionType& diag ) const
    {
      // extract diagonal entries 
      matrix().extractDiagonal( diag.blockVector() );
    }

    //! we only have right precondition
    bool rightPrecondition() const { return true; }

    //! precondition method for OEM Solvers 
    //! not fast but works, double is copied to block vector 
    //! and after application copied back
    void precondition(const double* arg, double* dest) const
    {
      createBlockVectors();
      
      assert( Arg_ );
      assert( Dest_ );

      RowBlockVectorType& Arg = *Arg_;
      ColumnBlockVectorType & Dest = *Dest_;
      
      // copy from double 
      double2Block(arg, Arg);
      
      // set Dest to zero 
      Dest = 0;
      
      assert( matrixAdap_ );
      // not parameter swaped for preconditioner 
      matrixAdap_->preconditionAdapter().apply(Dest , Arg);

      // copy back 
      block2Double( Dest , dest);
    }
      
    //! mult method for OEM Solver 
    void multOEM(const double* arg, double* dest) const
    {
      createBlockVectors();
      
      assert( Arg_ );
      assert( Dest_ );

      RowBlockVectorType& Arg = *Arg_;
      ColumnBlockVectorType & Dest = *Dest_;
      
      // copy from double 
      double2Block(arg, Arg);
      
      // call mult of matrix adapter  
      assert( matrixAdap_ );
      matrixAdap_->apply( Arg, Dest );

      //  copy back 
      block2Double( Dest , dest);
    }
    
    //! apply with discrete functions 
    void apply(const RowDiscreteFunctionType& arg,
               ColumnDiscreteFunctionType& dest) const 
    {
      createMatrixAdapter();
      assert( matrixAdap_ );
      matrixAdap_->apply( arg.blockVector(), dest.blockVector() );
    }

    //! apply with arbitrary discrete functions calls multOEM 
    template <class RowDFType, class ColDFType>
    void apply(const RowDFType& arg, ColDFType& dest) const 
    {
      multOEM( arg.leakPointer(), dest.leakPointer ());
    }

    //! mult method of matrix object used by oem solver
    void multOEM(const RowLeakPointerType& arg, ColumnLeakPointerType& dest) const
    {
      createMatrixAdapter();
      assert( matrixAdap_ );
      matrixAdap_->apply( arg.blockVector(), dest.blockVector() );
    }

    //! dot method for OEM Solver 
    double ddotOEM(const double* v, const double* w) const
    {
      createBlockVectors();
      
      assert( Arg_ );
      assert( Dest_ );

      RowBlockVectorType&    V = *Arg_;
      ColumnBlockVectorType& W = *Dest_;
      
      // copy from double 
      double2Block(v, V);
      double2Block(w, W);

#if HAVE_MPI 
      // in parallel use scalar product of discrete functions 
      ISTLBlockVectorDiscreteFunction< DomainSpaceType    > vF("ddotOEM:vF", domainSpace(), V ); 
      ISTLBlockVectorDiscreteFunction< RangeSpaceType > wF("ddotOEM:wF", rangeSpace(), W ); 
      return vF.scalarProductDofs( wF );
#else 
      return V * W;
#endif
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

    const DomainSpaceType& domainSpace() const { return domainSpace_; }
    const RangeSpaceType&  rangeSpace() const  { return rangeSpace_; }

    const RowMapperType& rowMapper() const { return rowMapper_; }
    const ColMapperType& colMapper() const { return colMapper_; }

    //! interface method from LocalMatrixFactory 
    ObjectType* newObject() const 
    {
      return new ObjectType(*this,
                            domainSpace(),
                            rangeSpace());
    }

    //! return local matrix object 
    LocalMatrixType localMatrix(const RowEntityType& rowEntity, 
                                const ColumnEntityType& colEntity) const 
    {
      return LocalMatrixType(localMatrixStack_,rowEntity,colEntity);
    }

  protected:  
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

    void removeObj( const bool alsoClearMatrix ) 
    {
      delete Dest_; Dest_ = 0;
      delete Arg_;  Arg_ = 0;
      delete matrixAdap_; matrixAdap_ = 0;
      
      if( alsoClearMatrix ) 
      {
        delete matrix_; 
        matrix_ = 0;
      }
    }

    // copy double to block vector 
    void double2Block(const double* arg, RowBlockVectorType& dest) const 
    {
      typedef typename RowBlockVectorType :: block_type BlockType;
      const size_t blocks = dest.size();
      int idx = 0;
      for(size_t i=0; i<blocks; ++i) 
      {
        BlockType& block = dest[i];
        enum { blockSize = BlockType :: dimension };
        for(int j=0; j<blockSize; ++j, ++idx) 
        {
          block[j] = arg[idx];
        }
      }
    }
      
    // copy block vector to double 
    void block2Double(const ColumnBlockVectorType& arg, double* dest) const 
    {
      typedef typename ColumnBlockVectorType :: block_type BlockType;
      const size_t blocks = arg.size();
      int idx = 0;
      for(size_t i=0; i<blocks; ++i) 
      {
        const BlockType& block = arg[i];
        enum { blockSize = BlockType :: dimension };
        for(int j=0; j<blockSize; ++j, ++idx) 
        {
          dest[idx] = block[j];
        }
      }
    }

    void createBlockVectors() const
    {
      if( ! Arg_ || ! Dest_ ) 
      { 
        delete Arg_; delete Dest_;
        Arg_  = new RowBlockVectorType( rowMapper_.size() );
        Dest_ = new ColumnBlockVectorType( colMapper_.size() );
      }

      createMatrixAdapter ();
    }

    void createMatrixAdapter () const 
    {
      if( ! matrixAdap_ ) 
      { 
        matrixAdap_ = new MatrixAdapterType(matrixAdapter());
      }
    }
    
  };



  // ISTLMatrixOperator
  // ------------------

  template< class DomainFunction, class RangeFunction, class TraitsImp >
  class ISTLMatrixOperator
  : public ISTLMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType, TraitsImp >,
    public Fem::AssembledOperator< DomainFunction, RangeFunction >
  {
    typedef ISTLMatrixOperator< DomainFunction, RangeFunction, TraitsImp > This;
    typedef ISTLMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType, TraitsImp > Base;

  public:
    typedef typename Base::DomainSpaceType DomainSpaceType;
    typedef typename Base::RangeSpaceType RangeSpaceType;

    /** \copydoc Fem::Operator::assembled */
    static const bool assembled = true ;

    using Base::apply;

    ISTLMatrixOperator ( const std::string &name,
                         const DomainSpaceType &domainSpace,
                         const RangeSpaceType &rangeSpace,
                         const std::string &paramfile = "" )
    : Base( domainSpace, rangeSpace, paramfile )
    {}

    virtual void operator() ( const DomainFunction &arg, RangeFunction &dest ) const
    {
      Base::apply( arg, dest );
    }

    const Base &systemMatrix () const
    {
      return *this;
    }
  };

  } // end namespace Fem

  // #if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: ISTLMatrixOperator ;
  using Fem :: ISTLMatrixObject ;

  // #endif // DUNE_FEM_COMPATIBILITY

} // end namespace Dune 

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_ISTLMATRIXWRAPPER_HH
