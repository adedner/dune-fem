namespace Dune
{

  namespace Fem
  {

    // where is the following used??
    // it is not used in this class => commented out
    // #define EPS 1.0E-15

    /*****************************/
    /*  Constructor(s)           */
    /*****************************/
    template <class T>
    SparseRowMatrix<T>::SparseRowMatrix(double omega) : omega_(omega)
    {
      values_ = 0;
      col_ = 0;
      dim_[0] = 0;
      dim_[1] = 0;
      memSize_ = 0;
      nz_ = 0;
      nonZeros_ = 0;
      if (checkNonConstMethods) assert(checkConsistency());
    }


    template <class T>
    SparseRowMatrix<T>::SparseRowMatrix(int rows, int cols, int nz,
                                        const T& dummy, double omega)
            : omega_(omega)
    {
      // standard settings as above
      values_ = 0;
      col_ = 0;
      dim_[0] = 0;
      dim_[1] = 0;
      memSize_ = 0;
      nz_ = 0;
      nonZeros_ = 0;

      // resize and get storage
      reserve(rows,cols,nz,dummy);

      // fill with value
      clear();

      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::removeObj()
    {
      if (checkNonConstMethods) assert(checkConsistency());
      if(values_) delete [] values_;
      if(col_) delete [] col_;
      if(nonZeros_) delete [] nonZeros_;
      values_ = 0;
      col_ = 0;
      nonZeros_ = 0;
      clearedRows_.clear();
      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    SparseRowMatrix<T>::~SparseRowMatrix()
    {
      if (checkNonConstMethods) assert(checkConsistency());
      removeObj();
      if (checkNonConstMethods) assert(checkConsistency());
    }

    /***********************************/
    /*  Construct from storage vectors */
    /***********************************/
    template <class T>
    void SparseRowMatrix<T>::
    reserve(int rows, int cols, int nz,const T& dummy )
    {
    // if (checkNonConstMethods) assert(checkConsistency());
      if( (rows == dim_[0]) && (cols == dim_[1]) && (nz == nz_))
      {
        clear();
        return;
      }

      removeObj();

      values_ = new T [ rows*nz ];
      col_    = new int [ rows*nz ];
      nonZeros_ = new int [ rows ];

      assert( values_ );
      assert( col_ );
      assert( nonZeros_ );

      dim_[0] = rows;
      dim_[1] = cols;

      memSize_ = rows * nz;
      nz_ = nz;
      // add first col for offset
      nz_ += firstCol ;

      assert( dim_[0] > 0 );
      assert( dim_[1] > 0 );

      // make resize
      newValues_.resize( nz_ );

      // only reserve for indices
      newIndices_.reserve( nz_ );

      // set all values to default values
      clear();
      if (checkNonConstMethods) assert(checkConsistency());

    }

    // resize with rows = cols = newSize
    template <class T>
    void SparseRowMatrix<T>::resize (int newSize)
    {
      resize(newSize,newSize);
    }

    // resize matrix
    template <class T>
    void SparseRowMatrix<T>::resize (int newRow, int newCol, int newNz )
    {
      if(newRow != this->size(0) || newNz > nz_ )
      {
        if( newNz < 0 ) newNz = nz_;

        int newMemSize = newRow * newNz ;

        int memHalf = (int) memSize_/2;
        if((newMemSize > memSize_) || (newMemSize < memHalf))
        {
          T tmp = 0;
          T * oldValues = values_;       values_ = 0;
          int * oldCol  = col_;          col_ = 0;
          int * oldNonZeros = nonZeros_; nonZeros_ = 0;
          const int oldNz = nz_;
          const int copySize = std::min( dim_[0] , newRow );
          const int oldSize = dim_[0];

          // reserve new memory
          reserve(newRow,newCol,newNz,tmp);

          if( (oldSize > 0) && (oldNz > 0 ))
          {
            std::memset( col_ , -1 , newRow * newNz * sizeof(int));
            for( int row = 0; row < copySize; ++ row )
            {
              const int newLoc = row * newNz ;
              const int oldLoc = row * oldNz ;
              std::memcpy( values_ + newLoc , oldValues + oldLoc , oldNz * sizeof(T) );
              std::memcpy( col_ + newLoc    , oldCol + oldLoc   , oldNz * sizeof(int) );
            }
            std::memcpy(nonZeros_, oldNonZeros, copySize * sizeof(int) );
          }

          delete [] oldValues;
          delete [] oldCol;
          delete [] oldNonZeros;
        }
        else
        {
          assert(newRow > 0);
          dim_[0] = newRow;
          dim_[1] = newCol;
        }
      }

      assert( this->size(0)  == newRow );
      assert( this->size(1)  == newCol );
    }

    //template <class T>
    //SparseRowMatrix<T>::
    //SparseRowMatrix(int rows, int cols, int nz, const T& val)
    //{
    //  reserve(rows,cols,nz,val);
    //}

    template< class T >
    inline T SparseRowMatrix<T>::operator() ( const int row, const int col ) const
    {
      assert( row >= 0 );
      assert( (row < dim_[0]) ? 1 : (std::cout << row << " bigger " << dim_[0] <<"\n", 0));

      const int nonZ = nonZeros_[row];
      int thisCol = row*nz_;
      for (int i=firstCol; i<nonZ; ++i)
      {
        if(col_[thisCol] == col)
        {
          return values_[thisCol];
        }
        ++thisCol;
      }
      return 0;
    }


    template< class T >
    inline T SparseRowMatrix<T>
      ::operator() ( const unsigned int row, const unsigned int col ) const
    {
      return (*this)( int( row ), int( col ) );
    }


    template <class T>
    int SparseRowMatrix<T>::colIndex(int row, int col)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      assert( row >= 0 );
      assert( row < dim_[0] );

      int i = 0;
      while ( i < nz_ && col_[row*nz_+i] < col && col_[row*nz_+i] != defaultCol )
        ++i;
      if (col_[row*nz_+i] == col)
        return i;  // column already in matrix
      else if ( col_[row*nz_+i] == defaultCol )
      { // add this column at end of this row
        ++nonZeros_[row];
        return i;
      }
      else
      {
        ++nonZeros_[row];
        // must shift this row to add col at the position i
        int j = nz_-1; // last column
        if (col_[row*nz_+j] != defaultCol)
        { // new space available - so resize
          resize( rows(), cols(), (2 * nz_) );
          j++;
        }
        for (;j>i;--j)
        {
          col_[row*nz_+j] = col_[row*nz_+j-1];
          values_[row*nz_+j] = values_[row*nz_+j-1];
        }
        col_[row*nz_+i] = col;
        values_[row*nz_+i] = 0;
        return i;
      }
    }

    template <class T>
    bool SparseRowMatrix<T>::find (int row, int col) const
    {
      int thisCol = 0;
      for(int i=firstCol; i<nz_; ++i)
      {
        thisCol = col_[row*nz_ +i];
        if(col == thisCol) return true;
        if(thisCol == defaultCol ) return false;
      }
      return false;
    }

    template <class T>
    void SparseRowMatrix<T>::clear()
    {
      T init = 0;
      for(register int i=0; i<dim_[0]*nz_; ++i)
      {
        values_ [i] = init;
        col_[i] = defaultCol;
      }

      for(register int i=0; i<dim_[0]; ++i)
      {
        nonZeros_[i] = 0;
      }

      clearedRows_.clear();

      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::clearRow(int row)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      assert( nonZeros_ );
      assert( values_ );
      assert( col_ );

      nonZeros_[row] = firstCol;

      int col = row * nz_;
      for(int i=0; i<nz_; ++i)
      {
        values_ [col] = 0;
        col_[col] = defaultCol;
        ++col;
      }

      // store row number (for UMFPACK solver)
      clearedRows_.insert( row );

      if (checkNonConstMethods) assert(checkConsistency());
    }


    template <class T>
    void SparseRowMatrix<T>::clearCol( int col )
    {
      if (checkNonConstMethods) assert(checkConsistency());
      assert( nonZeros_ );
      assert( values_ );
      assert( col_ );

      for(int i=0; i<dim_[0]; ++i)
        if((*this)(i,col)!= 0)
          set(i,col, 0);

      if (checkNonConstMethods) assert(checkConsistency());
    }



    template <class T>
    void SparseRowMatrix<T>::scaleRow(int row, const T& val )
    {
      assert( nonZeros_ );
      assert( values_ );
      assert( col_ );

      int col = row * nz_ ;
      for(int i=0; i<nz_ ; ++i, ++ col )
      {
        values_ [col] *= val ;
      }
    }

    template <class T>
    void SparseRowMatrix<T>::resort()
    {
      if (checkNonConstMethods) assert(checkConsistency());
      const int nRows = rows();
      for(int row=0; row<nRows; ++row)
      {
        resortRow(row);
      }
      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::resortRow(const int row)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      newIndices_.resize(0);
      int thisCol = row * nz_;

      for(int col=0; col<nz_; ++col)
      {
        int realCol =  col_[ thisCol + col ] ;
        if( realCol > defaultCol )
        {
          newIndices_.push_back( realCol );
        }
      }

      // set number of non zeros for row
      const int nZero = newIndices_.size();

      // nonZeros should be already at right size
      assert( nonZeros_[row] == nZero );
      //nonZeros_[row] = nZero;
      //std::cout << "found nz = " << nZero << "\n";

      // make values cache efficient
      std::sort( newIndices_.begin(), newIndices_.end() );
      for(int col=0; col<nZero; ++col)
      {
        int val = col_[ thisCol + col ];
        T value = values_[ thisCol + col ];
        for(int j=0; j<nZero; ++j)
        {
          if( newIndices_[j] == val )
            newValues_[j] = value;
        }
      }

      for(int col=0; col<nZero; ++col)
      {
        values_[ thisCol ] = newValues_[col];
        col_[ thisCol ] = newIndices_[col];
        ++thisCol;
      }
      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::set(int row, int col, T val)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      assert((col>=0) && (col <= dim_[1]));
      assert((row>=0) && (row <= dim_[0]));

      int whichCol = colIndex(row,col);
      assert( whichCol != defaultCol );

      {
        values_[row*nz_ + whichCol] = val;
        if (whichCol >= nonZeros_[row])
            nonZeros_[row]++;
        col_[row*nz_ + whichCol] = col;
      }
      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::add(int row, int col, T val)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      int whichCol = colIndex(row,col);
      assert( whichCol != defaultCol );
      values_[row*nz_ + whichCol] += val;
      col_[row*nz_ + whichCol] = col;
      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::multScalar(int row, int col, T val)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      int whichCol = colIndex(row,col);
      assert( whichCol != defaultCol );
      values_[row*nz_ + whichCol] *= val;
      col_[row*nz_ + whichCol] = col;
      if (checkNonConstMethods) assert(checkConsistency());
    }
    /***************************************/
    /*  Matrix-MV_Vector multiplication    */
    /***************************************/
    template <class T> template <class VECtype>
    void SparseRowMatrix<T>::mult(const VECtype *x, VECtype *ret) const
    {
      multOEM(x,ret);
    }

    template <class T> template <class VECtype>
    T SparseRowMatrix<T>::multOEMRow(const VECtype *x, const int row) const
    {
      T sum = 0;
      int thisCol = row*nz_ + firstCol ;
      const T * localValues = &values_[thisCol];
      const int nonZero = nonZeros_[row];
      for(int col = firstCol ; col<nonZero; ++col)
      {
        int realCol = col_[ thisCol ];
        assert( realCol > defaultCol );
        sum += localValues[col] * x[ realCol ];
        ++thisCol;
      }
      return sum;
    }

    template <class T> template <class VECtype>
    void SparseRowMatrix<T>::multOEM(const VECtype *x, VECtype *ret) const
    {
      for(register int row=0; row<dim_[0]; ++row)
      {
        ret[row] = multOEMRow( x, row );
      }
      return;
    }

    template <class T> template <class VECtype>
    void SparseRowMatrix<T>::multOEMAdd(const VECtype *x, VECtype *ret) const
    {
      for(register int row=0; row<dim_[0]; ++row)
      {
        ret[row] += multOEMRow( x, row );
      }
      return;
    }

    template <class T> template <class VECtype>
    void SparseRowMatrix<T>::multOEM_t(const VECtype *x, VECtype *ret) const
    {
      for(register int col=0; col<dim_[1]; ++col)
      {
        ret[col] = 0.0;
      }

      for(register int row=0; row<dim_[0]; ++row)
      {
        const int nonZero = nonZeros_[row];
        for(register int col=0; col<nonZero; ++col)
        {
          int thisCol = row*nz_ + col;
          int realCol = col_[ thisCol ];
          assert( realCol > defaultCol );
          ret[realCol] += values_[thisCol] * x[ row ];
        }
      }
      return;
    }

    /***************************************/
    /*  Matrix-MV_Vector multiplication    */
    /***************************************/

    template <class T> template <class ArgDFType, class DestDFType>
    void SparseRowMatrix<T>::apply(const ArgDFType &f, DestDFType &ret) const
    {
      typedef typename DestDFType::DofIteratorType DofIteratorType;

      typedef typename ArgDFType :: ConstDofBlockPtrType ConstDofBlockPtrType;
      enum { blockSize = ArgDFType :: DiscreteFunctionSpaceType :: localBlockSize };

      //! we assume that the dimension of the functionspace of f is the same as
      //! the size of the matrix
      DofIteratorType ret_it = ret.dbegin();

      for(int row=0; row<dim_[0]; ++row)
      {
        (*ret_it) = 0.0;

        //! DofIteratorType schould be the same
        for(int col=firstCol; col<nz_; ++col)
        {
          const int thisCol = row*nz_ + col;
          const int realCol = col_[thisCol];

          if( realCol == defaultCol ) continue;

          const int blockNr = realCol / blockSize ;
          const int dofNr = realCol % blockSize ;
          ConstDofBlockPtrType fBlock = f.block( blockNr );
          (*ret_it) += values_[thisCol] * (*fBlock)[ dofNr ];
        }

        ++ret_it;
      }
      return;
    }

    // apply to tranpose matrix
    template <class T> template <class ArgDFType, class DestDFType>
    void SparseRowMatrix<T>::apply_t(const ArgDFType &f, DestDFType &ret) const
    {
      typedef typename ArgDFType::ConstDofIteratorType ConstDofIteratorType;

      typedef typename DestDFType :: DofBlockPtrType DofBlockPtrType;
      enum { blockSize = DestDFType :: DiscreteFunctionSpaceType :: localBlockSize };

      //! we assume that the dimension of the functionspace of f is the same as
      //! the size of the matrix
      ret.clear();

      ConstDofIteratorType f_it = f.dbegin();

      for(int row=0; row<dim_[0]; ++row)
      {
        //! DofIteratorType schould be the same
        for(int col=firstCol; col<nz_; ++col)
        {
          const int thisCol = row * nz_ + col;
          const int realCol = col_[thisCol];

          if( realCol == defaultCol ) continue;

          const int blockNr = realCol / blockSize ;
          const int dofNr = realCol % blockSize ;
          DofBlockPtrType retBlock = ret.block( blockNr );

          (*retBlock)[ dofNr ] += values_[thisCol] * (*f_it);
        }

        ++f_it ;
      }
      return;
    }


    template <class T>
    void SparseRowMatrix<T>::print(std::ostream& s) const
    {
      s.precision( 6 );
      for(int row=0; row<dim_[0]; row++)
      {
        for(int col=0; col<dim_[1]; col++)
        {
          double val = (*this)(row,col);
          val = std::abs( val ) < 1e-14 ? 0 : val;
          s << val << " ";
        }
        s << "\n";
      }
      return;
    }

    template <class T>
    void SparseRowMatrix<T>::printReal(std::ostream& s) const
    {
      for(int row=0; row<dim_[0]; row++)
      {
        for(int col=0; col<nz_; col++)
        {
          s << values_[row*nz_ + col] << " ";
        }
        s << "\n";
      }
      return;
    }

    template <class T>
    void SparseRowMatrix<T>::printColumns(std::ostream& s) const
    {
      for(int row=0; row<dim_[0]; row++)
      {
        for(int col=0; col<nz_; col++)
        {
          s << col_[row*nz_ + col] << " ";
        }
        s << "\n";
      }
      return;
    }

    template <class T>
    void SparseRowMatrix<T>::printNonZeros(std::ostream& s) const
    {
      for(int row=0; row<dim_[0]; row++)
      {
        s << nonZeros_[row] << " ";
      }
      s << "\n";
      return;
    }

    template <class T>
    bool SparseRowMatrix<T>::checkConsistency() const
    {
      // check, whether nonzeros per row indeed correspond to reasonable
      // column-entries

      bool consistent = true;

      // only perform check, if there is any data:
      if (nonZeros_ || values_ || col_)
      {
        for(int row=0; row<dim_[0]; row++)
        {
          if (nonZeros_[row]<0 || nonZeros_[row]> dim_[1])
          {
            std::cout << "error in consistency of row " << row <<
                ": NonZeros_[row] = "<< nonZeros_[row] <<
                " is not within reasonable range "
                      << " of dim = (" << dim_[0]<<","<< dim_[1]<< ")\n";
            consistent = false;
            return(consistent);
          }

          for (int fakeCol =0; fakeCol < nonZeros_[row]; fakeCol++)
              if ((realCol(row,fakeCol)<0) || (realCol(row,fakeCol)>=dim_[1]))
              {

                std::cout << "error in consistency of row " << row <<
                    ": NonZeros_[row] = "<< nonZeros_[row] <<
                    ", fakeCol = " << fakeCol << ", realCol(fakeCol) = "
                          << realCol(row,fakeCol) << "\n" ;
                consistent = false;
                return(consistent);
              }
        }

        return(consistent);
      }

    //  std::cout << "consistent = " << consistent << "\n";

      assert(consistent);

      return(consistent);
    }

    template <class T>
    void SparseRowMatrix<T>::kroneckerKill(int row, int col)
    {
#ifndef NDEBUG
      if (checkNonConstMethods)
      {
        assert(checkConsistency());
      }
#endif
      unitRow(row);
      unitCol(col);
#ifndef NDEBUG
      if (checkNonConstMethods) assert(checkConsistency());
#endif
    }

    template <class T>
    void SparseRowMatrix<T>::unitRow(int row)
    {
      // only works for n x n matrices
      assert( dim_[ 0 ]  == dim_[ 1 ] );

      if (checkNonConstMethods) assert(checkConsistency());
      for(int i=1; i<nz_; i++)
      {
        values_[row*nz_ + i] = 0.0;
        col_[row*nz_ + i] = defaultCol;
      }
      values_[row*nz_] = 1.0;
      col_[row*nz_] = row;
      nonZeros_[row] = 1;

      // store row number (for UMFPACK solver)
      clearedRows_.insert( row );

      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::unitCol(int col)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      for(int i=0; i<dim_[0]; i++)
          if (i != col)
          {
            // only set 0 if nonzero column entry exists in current row
            if (this->operator()(i,col) != 0)
                set(i,col,0);
          }
          else set(col,col,1.0);
      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::checkSym()
    {
      if (checkNonConstMethods) assert(checkConsistency());
      double val;
      for(int i=0; i<this->size(0); i++)
      {
        for(int j=0; j<this->size(0); j++)
        {
          val = this->operator() (i,j);
          if(std::abs(val - this->operator() (j,i)) > 1E-10)
          {
            std::cout << val << " " << this->operator() (j,i) << " val \n";
          }
        }
      }
      if (checkNonConstMethods) assert(checkConsistency());
    }

    // diagonal conditioning
    template <class T> template <class DiscFuncType>
    void SparseRowMatrix<T>::getDiag(const ThisType & mass,
                                     const ThisType & B,
                                     DiscFuncType &diag) const
    {
      typedef typename DiscFuncType::DofIteratorType DofIteratorType;

      //! we assume that the dimension of the functionspace of f is the same as
      //! the size of the matrix
      DofIteratorType it = diag.dbegin();

      assert( this->size(0) == B.size(1) );
      assert( mass.size(0) == mass.size(1) );
      assert( mass.size(0) == this->size(1) );

      for(int row=0; row<this->size(0); ++row)
      {
        T sum = 0.0;
        int thisCol = row*nz_;
        const T * localValues = &values_[thisCol];
        const int nZeros = nonZeros_[row];
        for(register int col=0; col<nZeros; ++col)
        {
          int realCol = col_[ thisCol ];
          assert( realCol != defaultCol );
          double diag = mass(realCol,realCol);
          sum += diag * localValues[col] * B(realCol,row);
          ++thisCol;
        }

        (*it) = sum;
        ++it;
      }
      return;
    }

    // diagonal conditioning
    template <class T> template <class DiscFuncType>
    void SparseRowMatrix<T>::getDiag(const ThisType & B,
                                     DiscFuncType &diag) const
    {
      typedef typename DiscFuncType::DofIteratorType DofIteratorType;

      //! we assume that the dimension of the functionspace of f is the same as
      //! the size of the matrix
      DofIteratorType it = diag.dbegin();

      assert( this->size(0) == B.size(1) );

      for(int row=0; row<this->size(0); row++)
      {
        T sum = 0.0;
        for(register int col=0; col<nz_; ++col)
        {
          int thisCol = row*nz_ + col;
          int realCol = col_[ thisCol ];
          if ( realCol < 0 ) continue;
          sum += values_[ thisCol ] * B(realCol,row);
        }

        (*it) = sum;
        ++it;
      }
      return;
    }

    // diagonal conditioning
    template <class T> template <class DiscFuncType>
    void SparseRowMatrix<T>::getDiag(DiscFuncType &diag) const
    {
      typedef typename DiscFuncType::DofIteratorType DofIteratorType;

      //! we assume that the dimension of the functionspace of f is the same as
      //! the size of the matrix
      DofIteratorType it = diag.dbegin();

      for(int row=0; row<this->size(0); row++)
      {
        (*it) = (*this)(row,row);
        ++it;
      }
      return;
    }

    // add diagonal to given discrete function
    template <class T> template <class DiscFuncType>
    void SparseRowMatrix<T>::addDiag(DiscFuncType &diag) const
    {
      typedef typename DiscFuncType::DofIteratorType DofIteratorType;

      //! we assume that the dimension of the functionspace of f is the same as
      //! the size of the matrix
      DofIteratorType it = diag.dbegin();

      for(int row=0; row<this->size(0); row++)
      {
        (*it) += (*this)(row,row);
        ++it;
      }
      return;
    }

    template <class T>
    void SparseRowMatrix<T>::multiply(const SparseRowMatrix<T> & B,
        SparseRowMatrix<T> & res) const
    {
      //res.resize( B.size(0) );
      res.clear();
      //assert( res.numNonZeros() == B.numNonZeros() );

      //std::cout << res.numNonZeros() << "\n";
      for(int row=0; row<this->size(0); row++)
      {
        for(int col=0; col<B.size(1); col++)
        {
          T sum = 0;
          for(int k=0; k<B.size(0); k++)
          {
            sum += (*this)(row,k) * B(k,col);
          }

          if( std::abs(sum) > 0.0)
          {
            //std::cout << "add ("<<row<<","<<col<<")\n";
            res.add(row,col,sum);
          }
        }
      }

      //res.print(std::cout);
    }

    template <class T>
    void SparseRowMatrix<T>::add(const SparseRowMatrix<T> & B)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      assert( this->size(0) == B.size(0) );
      assert( nz_ == B.nz_ );

      for(register int i=0; i<dim_[0]*nz_; ++i)
      {
        assert( col_ [i] == B.col_ [i] );
        values_ [i] += B.values_[i];
      }
      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::scale(const T& factor)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      for(register int i=0; i<dim_[0]*nz_; ++i)
      {
        values_ [i] *= factor;
      }
      if (checkNonConstMethods) assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::ssorPrecondition(const T* u, T* x) const
    {
      const double omega = omega_;

      // (D - omega E) x = x_old (=u)
      for(int row=0; row<dim_[0]; ++row)
      {
        double diag=1.0, dot=0.0;
        // get row stuff
        int thisCol = row*nz_ + firstCol ;
        const T * localValues = &values_[thisCol];
        const int nonZero = nonZeros_[row];
        for(int col = firstCol ; col<nonZero; ++col)
        {
          const int realCol = col_[ thisCol ];
          assert( realCol > defaultCol );

          if (realCol < row)
          {
            dot += localValues[col] * x[realCol];
          }
          else if (realCol == row)
          {
            diag = localValues[col];
            assert( std::abs(diag) > 0.0 );
          }
          ++thisCol;
        }

        x[row] = (u[row] - omega*dot) / diag;
      }

      // D^{-1} (D - omega F) x = x_old (=x)
      for(int row=dim_[0]-1; row>=0; --row)
      {
        double diag=1.0, dot=0.0;
        int thisCol = row*nz_ + firstCol ;
        const T * localValues = &values_[thisCol];
        const int nonZero = nonZeros_[row];
        for(int col = firstCol ; col<nonZero; ++col)
        {
          const int realCol = col_[ thisCol ];
          assert( realCol > defaultCol );

          if (realCol > row)
          {
            dot += localValues[col] * x[realCol];
          }
          else if (realCol == row)
          {
            diag = localValues[col];
            assert( std::abs(diag) > 0.0 );
          }
          ++thisCol;
        }
        x[row] -= omega * dot / diag;
      }
    }

    template <class T>
    void SparseRowMatrix<T>::setupUMF(int n, int nAll, int* Ap, int* Ai, T* Ax,int &ANZ, int &LNZ)
    {
      // clear all columns that have been cleared from unitRow
      typedef typename std::set<int> :: iterator iterator ;
      const iterator end = clearedRows_.end();
      // needed for Lagrange dirichlet nodes
      if( clearedRows_.size() > 0 )
      {
        for(int row = 0; row < dim_[0]; ++ row)
        {
          if( clearedRows_.find( row ) != end )
          {
            unitCol( row );
          }
        }
      }

      int nZ = 0;
      unsigned int row = 0;
      for(int i=0; i<nAll; ++i)
      {
        if( i % nz_ == 0 )
        {
          const int row = (int)i / nz_ ;
          Ap[ row ] = nZ;
        }
        if( col_[i] != defaultCol && values_[i] != 0 )
        {
          Ai[ nZ ] = col_   [ i ];
          Ax[ nZ ] = values_[ i ];
          ++ nZ ;
          if (col_[i]<(int)row)
            ++ LNZ;
          else
            ++ ANZ;
        }
      }
      Ap[ n ] = nZ;
    }

    // use new method with symmetric/nonsymmetric flag
    template <class T>
    void SparseRowMatrix<T>::solveUMF(const T* b, T* x)
    {
#ifdef ENABLE_UMFPACK
#if 0 // LDL at the moment leads to unresolved symbols during linking
      int n = dim_[0];
      int nAll = n * nz_ ;
      int* Ap = new int [n+1];
      int* Ai = new int [ nAll ];
      T*   Ax = new   T [ nAll ];
      int ANZ,LNZ;
      setupUMF(n,nAll,Ap,Ai,Ax,ANZ,LNZ);

      double Lx [LNZ], D [n], Y [n] ;
      int Li [LNZ], Lp [n+1], Parent [n], Lnz [n], Flag [n], Pattern [n];
      /* factorize A into LDL’ (P and Pinv not used) */
      ldl_symbolic (n, Ap, Ai, Lp, Parent, Lnz, Flag, NULL, NULL) ;
      ldl_numeric (n, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D, Y, Pattern, Flag, NULL, NULL) ;
      /* solve Ax=b, overwriting b with the solution x so start with x=b*/
      for (int i=0;i<n;++i) x[i]=b[i];
      ldl_lsolve (n, x, Lp, Li, Lx) ;
      ldl_dsolve (n, x, D) ;
      ldl_ltsolve (n, x, Lp, Li, Lx) ;

      // delete temp memory
      delete [] Ap;
      delete [] Ax;
      delete [] Ai;
#else // want to use LDL but use umfpack at the moment
      const int n = dim_[0];
      const int nAll = n * nz_ ;

      int* Ap = new int [n+1];
      int* Ai = new int [ nAll ];
      T*   Ax = new   T [ nAll ];

      int ANZ,LNZ;
      setupUMF(n,nAll,Ap,Ai,Ax,ANZ,LNZ);

      int status;

      void *Symbolic, *Numeric;
      // symbolic analysis
      status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
      if (status != UMFPACK_OK)
      {
        if (status == UMFPACK_WARNING_singular_matrix)
          fprintf(stderr, "matrix is singluar!\n");
        else if(status == UMFPACK_ERROR_invalid_matrix)
          fprintf(stderr, "Number of entries in the matrix is negative, Ap [0] is nonzero, a column has a negative number of entries, a row index is out of bounds, or the columns of input matrix were jumbled (unsorted columns or duplicate entries).\n");
        else if(status == UMFPACK_ERROR_out_of_memory)
          fprintf(stderr, "Insufficient memory to perform the symbolic analysis.  If the analysis requires more than 2GB of memory and you are using the 32-bit (\"int\") version of UMFPACK, then you are guaranteed   to run out of memory.  Try using the 64-bit version of UMFPACK.\n");
        else if(status == UMFPACK_ERROR_argument_missing)
          fprintf(stderr, "One or more required arguments is missing.\n");
        else if(status == UMFPACK_ERROR_internal_error)
          fprintf(stderr, "omething very serious went wrong.  This is a bug.  Please contact the author (DrTimothyAldenDavis@gmail.com).\n");
        else
          fprintf(stderr, "umfpack_di_numeric() failed, %d\n", status);
      }
      // numeric analysis
      status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
      if (status != UMFPACK_OK) {
        if (status == UMFPACK_WARNING_singular_matrix)
          fprintf(stderr, "matrix is singluar!\n");
        else if(status == UMFPACK_ERROR_invalid_matrix)
          fprintf(stderr, "Number of entries in the matrix is negative, Ap [0] is nonzero, a column has a negative number of entries, a row index is out of bounds, or the columns of input matrix were jumbled (unsorted columns or duplicate entries).\n");
        else if(status == UMFPACK_ERROR_out_of_memory)
          fprintf(stderr, "Insufficient memory to perform the symbolic analysis.  If the analysis requires more than 2GB of memory and you are using the 32-bit (\"int\") version of UMFPACK, then you are guaranteed     to run out of memory.  Try using the 64-bit version of UMFPACK.\n");
        else if(status == UMFPACK_ERROR_argument_missing)
          fprintf(stderr, "One or more required arguments is missing.\n");
        else if(status == UMFPACK_ERROR_internal_error)
          fprintf(stderr, "omething very serious went wrong.  This is a bug. Please contact the author (DrTimothyAldenDavis@gmail.com).\n");
        else
          fprintf(stderr, "umfpack_di_numeric() failed, %d\n", status);
      }
      umfpack_di_free_symbolic (&Symbolic) ;

      // solve Ax = b
      // solve A^T x = b (since UMFPACK needs column wise storage, we got
      // row wise storage )
      status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
      if (status != UMFPACK_OK) {
        if (status == UMFPACK_WARNING_singular_matrix)
          fprintf(stderr, "matrix is singluar!\n");
        else if(status == UMFPACK_ERROR_invalid_matrix)
          fprintf(stderr, "Number of entries in the matrix is negative, Ap [0] is nonzero, a column has a negative number of entries, a row index is out of bounds, or the columns of input matrix were jumbled (unsorted columns or duplicate entries).\n");
        else if(status == UMFPACK_ERROR_out_of_memory)
          fprintf(stderr, "Insufficient memory to perform the symbolic analysis.  If the analysis requires more than 2GB of memory and you are using the 32-bit (\"int\") version of UMFPACK, then you are guaranteed     to run out of memory.  Try using the 64-bit version of UMFPACK.\n");
        else if(status == UMFPACK_ERROR_argument_missing)
          fprintf(stderr, "One or more required arguments is missing.\n");
        else if(status == UMFPACK_ERROR_internal_error)
          fprintf(stderr, "omething very serious went wrong.  This is a bug.  Please contact the author (DrTimothyAldenDavis@gmail.com).\n");
        else
          fprintf(stderr, "umfpack_di_numeric() failed, %d\n", status);
      }
      umfpack_di_free_numeric (&Numeric) ;

      // delete temp memory
      delete [] Ap;
      delete [] Ax;
      delete [] Ai;
#endif
#endif // ENABLE_UMFPACK
    }
    template <class T>
    void SparseRowMatrix<T>::solveUMFNonSymmetric(const T* b, T* x)
    {
#ifdef ENABLE_UMFPACK
      const int n = dim_[0];
      const int nAll = n * nz_ ;

      int* Ap = new int [n+1];
      int* Ai = new int [ nAll ];
      T*   Ax = new   T [ nAll ];

      int ANZ,LNZ;
      setupUMF(n,nAll,Ap,Ai,Ax,ANZ,LNZ);

      int status;

      // for nonsymmetric matrix we need to transpose A (column vs. row storagei)
      int* Cp = Ap;
      int* Ci = Ai;
      T* Cx = Ax;
      Ap = new int [n+1];
      Ai = new int [ nAll ];
      Ax = new   T [ nAll ];
      status = umfpack_di_transpose (n, n, Cp, Ci, Cx, (int *) NULL, (int *) NULL, Ap, Ai, Ax) ;
      delete [] Cx;
      delete [] Ci;
      delete [] Cp;

      void *Symbolic, *Numeric;
      // symbolic analysis
      status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);
      if (status != UMFPACK_OK)
      {
        if (status == UMFPACK_WARNING_singular_matrix)
          fprintf(stderr, "matrix is singluar!\n");
        else if(status == UMFPACK_ERROR_invalid_matrix)
          fprintf(stderr, "Number of entries in the matrix is negative, Ap [0] is nonzero, a column has a negative number of entries, a row index is out of bounds, or the columns of input matrix were jumbled (unsorted columns or duplicate entries).\n");
        else if(status == UMFPACK_ERROR_out_of_memory)
          fprintf(stderr, "Insufficient memory to perform the symbolic analysis.  If the analysis requires more than 2GB of memory and you are using the 32-bit (\"int\") version of UMFPACK, then you are guaranteed   to run out of memory.  Try using the 64-bit version of UMFPACK.\n");
        else if(status == UMFPACK_ERROR_argument_missing)
          fprintf(stderr, "One or more required arguments is missing.\n");
        else if(status == UMFPACK_ERROR_internal_error)
          fprintf(stderr, "omething very serious went wrong.  This is a bug.  Please contact the author (DrTimothyAldenDavis@gmail.com).\n");
        else
          fprintf(stderr, "umfpack_di_numeric() failed, %d\n", status);
      }
      // numeric analysis
      status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
      if (status != UMFPACK_OK) {
        if (status == UMFPACK_WARNING_singular_matrix)
          fprintf(stderr, "matrix is singluar!\n");
        else if(status == UMFPACK_ERROR_invalid_matrix)
          fprintf(stderr, "Number of entries in the matrix is negative, Ap [0] is nonzero, a column has a negative number of entries, a row index is out of bounds, or the columns of input matrix were jumbled (unsorted columns or duplicate entries).\n");
        else if(status == UMFPACK_ERROR_out_of_memory)
          fprintf(stderr, "Insufficient memory to perform the symbolic analysis.  If the analysis requires more than 2GB of memory and you are using the 32-bit (\"int\") version of UMFPACK, then you are guaranteed     to run out of memory.  Try using the 64-bit version of UMFPACK.\n");
        else if(status == UMFPACK_ERROR_argument_missing)
          fprintf(stderr, "One or more required arguments is missing.\n");
        else if(status == UMFPACK_ERROR_internal_error)
          fprintf(stderr, "omething very serious went wrong.  This is a bug. Please contact the author (DrTimothyAldenDavis@gmail.com).\n");
        else
          fprintf(stderr, "umfpack_di_numeric() failed, %d\n", status);
      }
      umfpack_di_free_symbolic (&Symbolic) ;

      // solve Ax = b
      // solve A^T x = b (since UMFPACK needs column wise storage, we got
      // row wise storage )
      status = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
      if (status != UMFPACK_OK) {
        if (status == UMFPACK_WARNING_singular_matrix)
          fprintf(stderr, "matrix is singluar!\n");
        else if(status == UMFPACK_ERROR_invalid_matrix)
          fprintf(stderr, "Number of entries in the matrix is negative, Ap [0] is nonzero, a column has a negative number of entries, a row index is out of bounds, or the columns of input matrix were jumbled (unsorted columns or duplicate entries).\n");
        else if(status == UMFPACK_ERROR_out_of_memory)
          fprintf(stderr, "Insufficient memory to perform the symbolic analysis.  If the analysis requires more than 2GB of memory and you are using the 32-bit (\"int\") version of UMFPACK, then you are guaranteed     to run out of memory.  Try using the 64-bit version of UMFPACK.\n");
        else if(status == UMFPACK_ERROR_argument_missing)
          fprintf(stderr, "One or more required arguments is missing.\n");
        else if(status == UMFPACK_ERROR_internal_error)
          fprintf(stderr, "omething very serious went wrong.  This is a bug.  Please contact the author (DrTimothyAldenDavis@gmail.com).\n");
        else
          fprintf(stderr, "umfpack_di_numeric() failed, %d\n", status);
      }
      umfpack_di_free_numeric (&Numeric) ;

      // delete temp memory
      delete [] Ap;
      delete [] Ax;
      delete [] Ai;
#endif
    }

  } // namespace Fem

} // namespace Dune
