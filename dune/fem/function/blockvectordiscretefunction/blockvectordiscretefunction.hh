// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH
#define DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH

#include <string>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/localfunction/standardlocalfunction.hh>
#include <dune/fem/function/localfunction/localfunctionwrapper.hh>

#include <dune/fem/storage/envelope.hh>
#include <dune/fem/misc/threadmanager.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh> 

#include <dune/grid/common/genericreferenceelements.hh>

namespace Dune {
namespace Fem {


  // forward declaration
  template< typename DiscreteFunctionSpace, typename BlockVector >
  class BlockVectorDiscreteFunction;

  /** \class IsBlockVectorDiscreteFunction
  *  \brief Tag for discrete functions using block vectors
  *
  *  A discrete function using block vectors for its dof storage and calculations should inherit from
  *  this struct. For example, Dune::Fem::MatrixOperator recognizes discrete functions with block vectors
  *  only as such if they inherit from this tag. If they do, their method .dofVector() is used (which is 
  *  a block vector). If they do not, this indicates that they don't use block vectors and thus provide 
  *  a .leakPointer() method -  which is used by Dune::Fem::MatrixOperator in this case.
  *
  */
  struct IsBlockVectorDiscreteFunction {};


  /** \class BlockVectorDiscreteFunctionTraits
  *  \brief Traits class for a BlockVectorDiscreteFunction
  *
  *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
  *  \tparam  BlockVector             implementation class of the block vector
  */
  template< typename DiscreteFunctionSpace, typename BlockVector >
  struct BlockVectorDiscreteFunctionTraits
  {
    typedef BlockVectorDiscreteFunctionTraits< DiscreteFunctionSpace, BlockVector > ThisType;
    typedef BlockVector                                                             DofVectorType;

    typedef DiscreteFunctionSpace                                               DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::DomainType                      DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType                       RangeType;
    typedef BlockVectorDiscreteFunction< DiscreteFunctionSpace, BlockVector >   DiscreteFunctionType;
    typedef Dune::StandardLocalFunctionFactory< ThisType >                      LocalFunctionFactoryType;
    typedef Dune::LocalFunctionStack< LocalFunctionFactoryType >                LocalFunctionStorageType;
    typedef typename LocalFunctionStorageType::LocalFunctionType                LocalFunctionType;
    typedef typename DofVectorType::IteratorType                                DofIteratorType;
    typedef typename DofVectorType::ConstIteratorType                           ConstDofIteratorType;
    typedef typename DofVectorType::DofBlockType                                DofBlockType;
    typedef typename DofVectorType::ConstDofBlockType                           ConstDofBlockType;
    typedef ::Dune::Envelope<DofBlockType>                                      DofBlockPtrType; 
    typedef ::Dune::Envelope<ConstDofBlockType>                                 ConstDofBlockPtrType;
    typedef typename DiscreteFunctionSpaceType::BlockMapperType                 MapperType;
    typedef typename DofVectorType::FieldType                                   DofType;
  };

  /** \class BlockVectorDiscreteFunctionTraits
  *  \brief A discrete function which uses block vectors for its dof storage and dof calculations
  *
  *  \tparam  DiscreteFunctionSpace   space the discrete function lives in
  *  \tparam  BlockVector             implementation class of the block vector
  */
  template< typename DiscreteFunctionSpace, typename BlockVector >
  class BlockVectorDiscreteFunction
  : public ::Dune::DiscreteFunctionInterface< BlockVectorDiscreteFunctionTraits< DiscreteFunctionSpace, BlockVector > >,
    public IsBlockVectorDiscreteFunction
  {

    /*
     I didn't implement these methods of DiscreteFunctionDefault (deliberately):
        void print ( std :: ostream &out ) const;
     */
    typedef BlockVectorDiscreteFunction< DiscreteFunctionSpace, BlockVector >   ThisType;
    typedef ParallelScalarProduct< ThisType >                                   ScalarProductType;

  public:
    // ==================== Types

    //! the traits of ThisType
    typedef BlockVectorDiscreteFunctionTraits< DiscreteFunctionSpace, BlockVector >     TraitsType;
    //! type for the discrete function space this function lives in
    typedef DiscreteFunctionSpace                                                       DiscreteFunctionSpaceType;
    //! type for the class which implements the block vector
    typedef BlockVector                                                                 BlockVectorType;
    //! type for the class which implements the block vector (which is the dof vector)
    typedef BlockVectorType                                                             DofVectorType;
    //! type for the mapper which maps local to global indices
    typedef typename TraitsType::MapperType                                             MapperType;
    //! type of the fields the dofs live in
    typedef typename TraitsType::DofType                                                DofType;
    //! pointer to a block of dofs
    typedef typename TraitsType::DofBlockPtrType                                        DofBlockPtrType;
    //! pointer to a block of dofs, const version
    typedef typename TraitsType::ConstDofBlockPtrType                                   ConstDofBlockPtrType;
    //! type of a block of dofs
    typedef typename TraitsType::DofBlockType                                           DofBlockType;
    //! iterator type for iterate over the dofs
    typedef typename TraitsType::DofIteratorType                                        DofIteratorType;
    //! iterator type for iterate over the dofs, const verision
    typedef typename TraitsType::ConstDofIteratorType                                   ConstDofIteratorType;
    //! type of the range field
    typedef typename DiscreteFunctionSpaceType::RangeFieldType                          RangeFieldType;
    //! type of the domain field
    typedef typename DiscreteFunctionSpaceType::DomainFieldType                         DomainFieldType;
    //! type of the local functions
    typedef typename TraitsType::LocalFunctionType                                      LocalFunctionType;
    //! type of the local function storage
    typedef typename TraitsType::LocalFunctionStorageType                               LocalFunctionStorageType;
    //! type of the factory which produces local functions
    typedef typename TraitsType::LocalFunctionFactoryType                               LocalFunctionFactoryType;
    //! type of the discrete functions's domain
    typedef typename TraitsType::DomainType                                             DomainType;
    //! type of the discrete functions's range
    typedef typename TraitsType::RangeType                                              RangeType;
    //! size type of the block vector
    typedef typename BlockVectorType::SizeType                                          SizeType;

    //! size of the dof blocks
    enum { blockSize = BlockVectorType::blockSize };

    /** \brief Constructor to use if the vector storing the dofs (which is a block vector) already exists
     *
     *  \param[in]  name         name of the discrete function
     *  \param[in]  dfSpace      space the discrete function lives in
     *  \param[in]  blockVector  reference to the blockVector
     */
    BlockVectorDiscreteFunction ( const std::string& name,
                                  const DiscreteFunctionSpaceType& dfSpace,
                                  BlockVectorType& blockVector )
    : dfSpace_( dfSpace ),
      lfFactory_( *this ),
      #ifdef USE_SMP_PARALLEL
        lfStorageVec_ ( Fem :: ThreadManager :: maxThreads() ),
      #else 
        lfStorage_( lfFactory_ ),
      #endif
      mapper_( dfSpace.blockMapper() ),
      name_( name ),
      memPair_( static_cast< DofStorageInterface* >( 0 ), &blockVector ),
      scalarProduct_( dfSpace_ )
    {
      #ifdef USE_SMP_PARALLEL
        for( size_t i=0 ; i<lfStorageVec_.size(); ++i ) 
        {
          lfStorageVec_[ i ] = new LocalFunctionStorageType( lfFactory_ );
        }
      #endif
    }

    /** \brief Constructor to use if the vector storing the dofs does not exist yet
     *
     *  \param[in]  name         name of the discrete function
     *  \param[in]  dfSpace      space the discrete function lives in
     */
    BlockVectorDiscreteFunction ( const std::string &name,
                                  const DiscreteFunctionSpaceType &dfSpace )
    : dfSpace_( dfSpace ),
      lfFactory_( *this ), 
      #ifdef USE_SMP_PARALLEL
        lfStorageVec_ ( Fem :: ThreadManager :: maxThreads() ),
      #else 
        lfStorage_( lfFactory_ ),
      #endif
      mapper_( dfSpace.blockMapper() ),
      name_( name ),
      memPair_( allocateManagedDofStorage< BlockVectorType >( space().grid(), mapper_, name ) ),
      scalarProduct_( dfSpace_ )
    {
      #ifdef USE_SMP_PARALLEL
        for( size_t i=0 ; i<lfStorageVec_.size(); ++i ) 
        {
          lfStorageVec_[ i ] = new LocalFunctionStorageType( lfFactory_ );
        }
      #endif
    }


    // TODO: DiscreteFunctionDefault prohibits the copy ctor. Should this be done here too?

    /** \brief Copy constructor
     */
    BlockVectorDiscreteFunction ( const ThisType &other )
    : dfSpace_( other.space() ),
      lfFactory_( *this ),
      #ifdef USE_SMP_PARALLEL
        lfStorageVec_ ( Fem :: ThreadManager :: maxThreads() ),
      #else 
        lfStorage_( lfFactory_ ),
      #endif
      mapper_( other.space().blockMapper() ), 
      name_( other.name() ),
      memPair_( allocateManagedDofStorage< BlockVectorType >( space().grid(), mapper_, name() ) ),
      scalarProduct_( dfSpace_ )
    {
      dofVector() = other.dofVector();

      // we create our own factories
      #ifdef USE_SMP_PARALLEL
        for( size_t i=0 ; i<lfStorageVec_.size(); ++i ) 
        {
          lfStorageVec_[ i ] = new LocalFunctionStorageType( lfFactory_ );
        }
      #endif

    }

    ~BlockVectorDiscreteFunction ()
    {
      // TODO: use a smart pointer for this?
      // No need for a null check here. Stroustrup: "Applying delete to zero has no effect."
      delete memPair_.first;

      #ifdef USE_SMP_PARALLEL
        for( size_t i=0 ; i<lfStorageVec_.size(); ++i )
        {
          delete lfStorageVec_[ i ]; 
        }
      #endif
    }

  private:

    // an empty constructor would not make sense for a discrete function
    BlockVectorDiscreteFunction ();

    // TODO: un-disallow this??
    ThisType& operator= (const ThisType& other);

  public:

    /** \brief Copy other's dof vector to *this
     *
     *  \param[in]  other   reference to the other dof vector
     *  \return Reference to this
     */
    void assign (const ThisType& other)
    {
      dofVector().assign(other.dofVector());
    }
  
    /** \brief Add scalar*v to *this
     *
     *  \param[in]  scalar  scalar by which v has to be multiplied before adding it to *this
     *  \param[in]  v       the other discrete function which has to be scaled and added
     */
    void axpy (const RangeFieldType& scalar, const ThisType& v)
    {
      dofVector().addScaled(v.dofVector(), scalar);
    }

    /** \brief Add scalar*v to *this
     *
     *  \param[in]  scalar  scalar by which v has to be multiplied before adding it to *this
     *  \param[in]  v       the other discrete function which has to be scaled and added
     */
    void axpy (const ThisType& v, const RangeFieldType& scalar)
    {
      dofVector().addScaled(v.dofVector(), scalar);
    }

    /** \brief This method does the same as axpy
     *   These two methods are here to make this class compatible to the old interface - the newton solver relies on this method for example...
     *   In new code, axpy should be used.
     *
     *  \param[in]  scalar  scalar by which v has to be multiplied before adding it to *this
     *  \param[in]  v       the other discrete function which has to be scaled and added
     */
    void addScaled (const RangeFieldType& scalar, const ThisType& v) { axpy(scalar, v); } // TODO: make these methods DUNE_DEPRECATED?
    void addScaled (const ThisType& v, const RangeFieldType& scalar) { axpy(scalar, v); }

    /** \brief Obtain the (modifiable) 'index'-th block
     *
     *  \param[in]  index   index of the block
     *  \return The (modifiable) 'index'-th block
     */
    DofBlockPtrType block ( unsigned int index )
    {
      return DofBlockPtrType( memPair_.second->operator[](index) );
    }

    /** \brief Obtain the (constant) 'index'-th block
     *
     *  \param[in]  index   index of the block
     *  \return The (constant) 'index'-th block
     */
    ConstDofBlockPtrType block ( unsigned int index ) const 
    {
      return ConstDofBlockPtrType( memPair_.second->operator[](index) );
    }

    /** \brief Set each dof to zero
     */
    void clear ()
    {
       dofVector().clear();
    }

    /** \copydoc Dune::DiscreteFunctionInterface::communicate() */
    void communicate()
    {
      // Works only in singleThreadMode currently
      assert( Fem :: ThreadManager :: singleThreadMode() );
      space().communicate( *this );
    }
  
    /** \copydoc Dune::DiscreteFunctionInterface::dataHandle */
    template< class Operation >
    typename DiscreteFunctionSpaceType :: template CommDataHandle< ThisType, Operation > :: Type 
    dataHandle ( const Operation *operation )
    {
      return space().createDataHandle( *this, operation );
    }

    /** \brief Obtain the constant iterator pointing to the first dof
     *
     *  \return Constant iterator pointing to the first dof
     */
    ConstDofIteratorType dbegin () const { return dofVector().dbegin(); }

    /** \brief Obtain the non-constant iterator pointing to the first dof
     *
     *  \return Non-Constant iterator pointing to the first dof
     */
    DofIteratorType dbegin () { return dofVector().dbegin(); }

    /** \brief Obtain the constant iterator pointing to the last dof
     *
     *  \return Constant iterator pointing to the last dof
     */
    ConstDofIteratorType dend () const { return dofVector().dend(); }

    /** \brief Obtain the non-constant iterator pointing to the last dof
     *
     *  \return Non-Constant iterator pointing to the last dof
     */
    DofIteratorType dend () { return dofVector().dend(); }


    /** \copydoc Dune::DiscreteFunctionInterface::dofsValid */
    bool dofsValid () const
    {
      // check for NaN or inf...
      const ConstDofIteratorType end = dend();
      for( ConstDofIteratorType it = dbegin(); it != end; ++it )
        if( *it != *it )
          return false;

      return true;
    }

    /** \brief Obtain constant reference to the dof vector
     *
     *  \return Constant reference to the block vector
     */
    const BlockVectorType &dofVector () const
    {
      assert( memPair_.second );
      return *memPair_.second;
    }

    /** \brief Obtain reference to the dof vector
     *
     *  \return Reference to the block vector
     */
    BlockVectorType &dofVector ()
    {
      assert( memPair_.second );
      return *memPair_.second;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::enableDofCompression()
     */
    void enableDofCompression ()
    {
      if ( memPair_.first )
        memPair_.first->enableDofCompression();
    }

    /** \copydoc Dune::Fem::Function::evaluate(const DomainType &x,RangeType &ret) const */
    void evaluate ( const DomainType &x, RangeType &ret ) const
    {
      FieldVector< int, 0 > diffVariable;
      evaluate( diffVariable, x, ret );
    }

    /** \copydoc Dune::Fem::Function::evaluate(const FieldVector<int,diffOrder> &diffVariable,const DomainType &x,RangeType &ret) const */
    template< int diffOrder >
    void evaluate ( const FieldVector< int, diffOrder > &diffVariable,
                    const DomainType &x,
                    RangeType &ret ) const
    {
      //std::cout << "using my evaluate() :)\n\n\n\n";
      typedef typename DiscreteFunctionSpaceType::IteratorType Iterator;
      typedef typename Iterator::Entity Entity;
      typedef typename Entity::Geometry Geometry;
      typedef typename Geometry :: LocalCoordinate LocalCoordinateType;

      const int dimLocal = LocalCoordinateType :: dimension;
      
      const Iterator end = space().end();
      for( Iterator it = space().begin(); it != end; ++it )
      {
        const Entity &entity = *it;
        const Geometry &geometry = entity.geometry();

        const GenericReferenceElement< DomainFieldType, dimLocal > &refElement
          = GenericReferenceElements< DomainFieldType, dimLocal >::general( geometry.type() );

        const LocalCoordinateType xlocal = geometry.local( x );
        if( refElement.checkInside( xlocal ) )
        {
          localFunction( entity ).evaluate( diffVariable, xlocal, ret );
          return;
        }
      }
      DUNE_THROW( RangeError, "DiscreteFunctionDefault::evaluate: x is not within domain." );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) */
    template< class EntityType >
    LocalFunctionType localFunction ( const EntityType &entity )
    {
      return localFunctionStorage().localFunction( entity );
    }

    /** \copydoc Dune::DiscreteFunctionInterface::localFunction(const EntityType &entity) const */
    template< class EntityType >
    const LocalFunctionType localFunction ( const EntityType &entity ) const
    {
      return localFunctionStorage().localFunction( entity );
    }


    /** \copydoc Dune::DiscreteFunctionInterface::name() const */
    const std::string& name () const { return name_; }

    /** \brief Add another discrete function to this one
     *
     *  \param[in]  other   discrete function to add
     *  \return  constant reference to *this
     */
    const ThisType &operator+= ( const ThisType &other )
    {
      dofVector() += other.dofVector();
      return *this;
    }

    /** \brief Subtract another discrete function from this one
     *
     *  \param[in]  other   Discrete function to subtract
     *  \return Constand reference to this
     */
    const ThisType &operator-= ( const ThisType &other )
    {
      dofVector() -= other.dofVector();
      return *this;
    }

    /** \brief Scale this
     *
     *  \param[in] scalar   scalar factor for the scaling
     *  \return Constant reference to *this
     */
    const ThisType &operator*= ( const DofType &scalar )
    {
      dofVector() *= scalar;
      return *this;
    }

    /** \brief Divide each dof by a scalar
     *
     *  \param[in] scalar   Scalar to divide each dof by
     *  \return Constant reference to *this
     */
    const ThisType &operator/= ( const DofType &scalar )
    {
      dofVector() *= 1./scalar;
      return *this;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::read */
    template< class StreamTraits >
    void read ( InStreamInterface< StreamTraits > &in )
    {
      unsigned int versionId = in.readUnsignedInt();
      if( versionId < DUNE_VERSION_ID(0,9,1) )
        DUNE_THROW( IOError, "Trying to read outdated file." );
      else if( versionId > DUNE_MODULE_VERSION_ID(DUNE_FEM) )
        std :: cerr << "Warning: Reading discrete function from newer version: "
                    << versionId << std :: endl;

      in >> name_;
      
      if( static_cast< SizeType >( in.readInt() ) != size() && 
          size() != static_cast< SizeType >( space().size() ) ) // only read compressed vectors 
        DUNE_THROW( IOError, "Trying to read discrete function of different size." );

      const DofIteratorType end = dend();
      for( DofIteratorType it = dbegin(); it != end; ++it )
        in >> *it;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::read_xdr */
    bool read_xdr ( const std :: string filename )
    {
      // TODO: do we want to implement this?
      try
      {
        XDRFileInStream in( filename );
        read( in );
        return true;
      }
      catch( Exception e )
      {
        return false;
      }
    }

    /** \copydoc Dune::DiscreteFunctionInterface::read_ascii */
    bool read_ascii ( const std :: string filename )
    {
      // TODO: do we want to implement this?
      try
      {
        ASCIIInStream in( filename );
        read( in );
        return true;
      }
      catch( Exception e )
      {
        return false;
      }
    }
 
  
    /** \copydoc Dune::DiscreteFunctionInterface::scalarProductDofs */
    DofType scalarProductDofs ( const ThisType &other ) const
    {
      return scalarProduct_.scalarProductDofs( *this, other );
    }

    /** \brief Return the number of blocks in the block vector
     *
     *  \return Number of block in the block vector
     */
    SizeType size () const { return dofVector().size(); }

    /** \copydoc Dune::DiscreteFunctionInterface::space() const */
    const DiscreteFunctionSpaceType& space () const { return dfSpace_; }

    /** \copydoc Dune::DiscreteFunctionInterface::write */
    template< class StreamTraits >
    void write ( OutStreamInterface< StreamTraits > &out ) const
    {
      out << DUNE_MODULE_VERSION_ID(DUNE_FEM);
      out << name_;
    
      // only allow write when vector is compressed 
      if( size() != static_cast< SizeType >( space().size() ) )
        DUNE_THROW(InvalidStateException,"Writing DiscreteFunction in uncompressed state!");
      
      out << size();

      const ConstDofIteratorType end = dend();
      for( ConstDofIteratorType it = dbegin(); it != end; ++it )
        out << *it;
    }

    /** \copydoc Dune::DiscreteFunctionInterface::write_ascii */
    bool write_ascii ( const std :: string filename ) const
    {
      // TODO: do we want to implement this?
      try
      {
        ASCIIOutStream out( filename );
        write( out );
        return true;
      }
      catch( Exception e )
      {
        return false;
      }
    }

    /** \copydoc Dune::DiscreteFunctionInterface::write_xdr */
    bool write_xdr ( const std :: string filename ) const
    {
      // TODO: do we want to implement this here?
      try
      {
        XDRFileOutStream out( filename );
        write( out );
        return true;
      }
      catch( Exception e )
      {
        return false;
      }
    }

  private:

    /*
     * ============================== private methods ====================
     */

    // Obtain the correct local function storage object for this thread
    LocalFunctionStorageType& localFunctionStorage () const
    {
      #ifdef USE_SMP_PARALLEL
        return *(lfStorageVec_[ Fem :: ThreadManager :: thread() ]);
      #else 
        return lfStorage_;
      #endif
    }

    /* 
     * ============================== data fields ====================
     */
    const DiscreteFunctionSpaceType& dfSpace_;
    const LocalFunctionFactoryType lfFactory_;
    // local function storage
    #ifdef USE_SMP_PARALLEL
      mutable std::vector< LocalFunctionStorageType* > lfStorageVec_;
    #else 
      mutable LocalFunctionStorageType lfStorage_;
    #endif
    MapperType mapper_;
    std::string name_;
    std::pair< DofStorageInterface *, BlockVectorType * > memPair_;
    ScalarProductType scalarProduct_;

  };

} // namespace Fem
} // namespace Dune

#endif // #ifndef DUNE_FEM_BLOCKVECTORDISCRETEFUNCTION_HH
