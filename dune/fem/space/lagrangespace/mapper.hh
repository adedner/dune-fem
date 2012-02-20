#ifndef DUNE_LAGRANGESPACE_MAPPER_HH
#define DUNE_LAGRANGESPACE_MAPPER_HH

#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/fem/misc/capabilities.hh>
#include <dune/fem/misc/codimmap.hh>
#include <dune/fem/misc/metaprogramming.hh>
#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/dofmapper/localkey.hh>
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>
#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>

namespace Dune
{
  
  template< class GridPart, int polOrder >
  class LagrangeMapper;


  
  template< class GridPart, int polOrder >
  struct LagrangeMapperTraits
  {
    typedef GridPart GridPartType;
    
    static const int polynomialOrder = polOrder;

    typedef typename GridPartType::template Codim< 0 >::IteratorType::Entity EntityType;
    typedef LagrangeMapper< GridPartType, polynomialOrder > DofMapperType;
    typedef DefaultDofMapIterator< EntityType, DofMapperType > DofMapIteratorType;
  };



  // First Order Lagrange Mapper
  // ---------------------------

  template< class GridPart >
  class LagrangeMapper< GridPart, 1 >
  : public CodimensionMapper< GridPart, GridPart::GridType::dimension >
  {
    typedef LagrangeMapper< GridPart, 1 > ThisType;
    typedef CodimensionMapper< GridPart, GridPart::GridType::dimension > BaseType;

  public:
    //! type of the grid part
    typedef typename BaseType::GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename BaseType::EntityType EntityType;

    //! type of the underlying grid
    typedef typename GridPartType::GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType::ctype FieldType;

    //! dimension of the grid
    static const int dimension = GridType::dimension;

    //! order of the Lagrange polynoms
    static const int polynomialOrder = 1; 

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder > LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std::vector< const LagrangePointSetType * > LagrangePointSetContainerType;

  public:
    //! constructor
    LagrangeMapper ( const GridPartType &gridPart, const LagrangePointSetContainerType & )
    : BaseType( gridPart )
    {}

    bool fixedDataSize ( const int codim ) const
    {
      return true;
    }
  };



  // Second Order Lagrange Mapper
  // ----------------------------

  template< class GridPart >
  class LagrangeMapper< GridPart, 2 >
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, 2 > >
  {
    typedef LagrangeMapper< GridPart, 2 > ThisType;
    typedef DofMapperDefault< LagrangeMapperTraits< GridPart, 2 > > BaseType;

  public:
    typedef LagrangeMapperTraits< GridPart, 2 > Traits;
    
    //! type of the grid part
    typedef typename Traits::GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename Traits::EntityType EntityType;

    //! type of DofMapIterator
    typedef typename Traits::DofMapIteratorType DofMapIteratorType;
 
    //! type of the underlying grid
    typedef typename GridPartType::GridType GridType;

    //! type of the index set
    typedef typename GridPartType::IndexSetType IndexSetType;

    //! type of coordinates within the grid
    typedef typename GridType::ctype FieldType;

    //! dimension of the grid
    static const int dimension = GridType::dimension;

    //! order of the Lagrange polynoms
    static const int polynomialOrder = Traits::polynomialOrder;

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std::vector< const LagrangePointSetType * > LagrangePointSetContainerType;

    //! type of the DoF manager
    typedef DofManager< GridType > DofManagerType;

  public:
    //! constructor
    LagrangeMapper ( const GridPartType &gridPart,
                     const LagrangePointSetContainerType &lagrangePointSetContainer )
    : dm_( DofManagerType :: instance(gridPart.grid()) ),
      indexSet_( gridPart.indexSet() ),
      lagrangePointSetContainer_( lagrangePointSetContainer ),
      overShoot_( Parameter::getValidValue( "fem.lagrangemapper.overshoot", double( 1.5 ), ValidateNotLess< double >( 1.0 ) ) )
    {
      numDofs_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        maxDofs_[ codim ] = 0;
      
      typedef typename LagrangePointSetContainerType::const_iterator IteratorType;
      const IteratorType end = lagrangePointSetContainer_.end();
      for( IteratorType it = lagrangePointSetContainer_.begin(); it != end; ++it )
      {
        const LagrangePointSetType *lagrangePointSet = *it;
        if( lagrangePointSet )
        {
          numDofs_ = std::max( numDofs_, lagrangePointSet->size() );
          for( int codim = 0; codim <= dimension; ++codim )
            maxDofs_[ codim ] = std::max( maxDofs_[ codim ], lagrangePointSet->maxDofs( codim ) );
        }
      }

      computeOffsets();
      dm_.addIndexSet( *this );
    }
    
    //! destructor 
    virtual ~LagrangeMapper ()
    {
      dm_.removeIndexSet( *this );
    }

    //! return overall number of degrees of freedom 
    int size () const
    {
      return size_;
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      const Fem::LocalKey &dofInfo = lagrangePointSet( entity.type() ).localKey( localDof );
      const unsigned int codim = dofInfo.codim();
      return offset_[ codim ] + indexSet_.subIndex( entity, dofInfo.subEntity(), codim );
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      static const unsigned int codim = Entity::codimension;
      assert( (localDof >= 0) && (localDof < numEntityDofs( entity )) );
      return offset_[ codim ] + indexSet_.index( entity );
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return numDofs_;
    }

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      return lagrangePointSet( entity.type() ).size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return maxDofs_[ Entity::codimension ];
    }
    
    /** \brief Check, whether any DoFs are associated with a codimension */
    bool contains ( int codim ) const
    {
      return true;
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    bool fixedDataSize ( int codim ) const
    {
      return false;
    }

    /** \copydoc Dune::DofMapper::oldIndex */
    int oldIndex ( const int hole, const int codim ) const
    {
      return offset_[ codim ] + indexSet_.oldIndex( hole, codim );
    }

    /** \copydoc Dune::DofMapper::newIndex */
    int newIndex ( const int hole, const int codim ) const
    {
      return offset_[ codim ] + indexSet_.newIndex( hole, codim );
    }

    /** \copydoc Dune::DofMapper::numberOfHoles */
    int numberOfHoles ( const int codim ) const
    {
      return maxDofs_[ codim ] * indexSet_.numberOfHoles( codim );
    }

    /** \copydoc Dune::DofMapper::numBlocks
     */
    int numBlocks () const 
    {
      return dimension + 1;
    }

    /** \copydoc Dune::DofMapper::oldOffset
     */
    int oldOffSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return oldOffSet_[ block ];
    }

    /** \copydoc Dune::DofMapper::newOffset
     */
    int offSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return offset_[ block ];
    }

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return BaseType::checkConsecutive( indexSet_ );
    }


    // Adaptation Methods (as for Index Sets)

    void insertEntity ( const EntityType &entity )
    {
      // check, whether we need to enlarge any block and compute oversized block size
      unsigned int oldUpperBound = size_;
      for( int codim = dimension ; codim >= 0; --codim )
      {
        const unsigned int codimSize = indexSet_.size( codim ) * maxDofs_[ codim ];
        const unsigned int newUpperBound = offset_[ codim ] + codimSize;
        if( newUpperBound > oldUpperBound )
          return computeOffsets( overShoot_ );
        oldUpperBound = offset_[ codim ];
      }
    }

    void removeEntity ( const EntityType &entity )
    {}

    void resize ()
    {
      computeOffsets();
    }

    bool compress ()
    {
      computeOffsets();
      return true;
    }

    void read_xdr ( const char *filename, int timestep )
    {
      computeOffsets();
    }

    void write_xdr ( const char *filename, int timestep )
    {}

  private:
    // prohibit copying and assignment
    LagrangeMapper ( const ThisType & );
    ThisType &operator=( const ThisType & );

    void computeOffsets ( const double overShoot = 1.0 )
    {
      // store old offset and calculate new offsets
      size_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        oldOffSet_[ codim ] = offset_[ codim ];
        offset_[ codim ] = size_;

        const unsigned int codimSize = indexSet_.size( codim ) * maxDofs_[ codim ];
        size_ += static_cast< unsigned int >( overShoot * codimSize );
      }
    }

    const LagrangePointSetType &lagrangePointSet ( const GeometryType type ) const
    {
      const LagrangePointSetType *lagrangePointSet
        = lagrangePointSetContainer_[ LocalGeometryTypeIndex::index( type ) ];
      assert( lagrangePointSet );
      return *lagrangePointSet;
    }

    // reference to dof manager
    DofManagerType &dm_;

    const IndexSetType &indexSet_;
    
    const LagrangePointSetContainerType &lagrangePointSetContainer_;

    // memory overshoot 
    const double overShoot_ ;

    unsigned int maxDofs_[ dimension+1 ];
    mutable unsigned int offset_[ dimension+1 ];
    mutable unsigned int oldOffSet_[ dimension+1 ];
    mutable unsigned int size_;
    unsigned int numDofs_;
  };



  // Higher Order Lagrange Mapper
  // ----------------------------
  //
  // Note: This mapper assumes that the grid is "twist-free".

#ifdef USE_TWISTFREE_MAPPER
  template< class GridPart, int polOrder >
  class LagrangeMapper
  : public DofMapperDefault< LagrangeMapperTraits< GridPart, polOrder > >
  {
    typedef LagrangeMapper< GridPart, polOrder > ThisType;
    typedef DofMapperDefault< LagrangeMapperTraits< GridPart, polOrder > > BaseType;

  public:
    typedef LagrangeMapperTraits< GridPart, polOrder > Traits;
    
    //! type of the grid part
    typedef typename Traits::GridPartType GridPartType;

    //! type of entities (codim 0)
    typedef typename Traits::EntityType EntityType;

    //! type of DofMapIterator
    typedef typename Traits::DofMapIteratorType DofMapIteratorType;
 
    //! type of the underlying grid
    typedef typename GridPartType::GridType GridType;

    //! type of coordinates within the grid
    typedef typename GridType::ctype FieldType;

    //! dimension of the grid
    static const int dimension = GridType::dimension;

    //! order of the Lagrange polynoms
    static const int polynomialOrder = Traits::polynomialOrder;

    //! type of the index set
    typedef typename GridPartType::IndexSetType IndexSetType;

    //! type of the Lagrange point set
    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    //! type of the map for the Lagrange point sets
    typedef std::vector< const LagrangePointSetType * > LagrangePointSetContainerType;

    //! type of the DoF manager
    typedef DofManager< GridType > DofManagerType;

  public:
    //! constructor
    LagrangeMapper ( const GridPartType &gridPart,
                     const LagrangePointSetContainerType &lagrangePointSetContainer )
    : dm_( DofManagerType::instance(gridPart.grid()) ),
      indexSet_( gridPart.indexSet() ),
      lagrangePointSetContainer_( lagrangePointSetContainer ),
      overShoot_( Parameter::getValidValue( "fem.lagrangemapper.overshoot", double( 1.5 ), ValidateNotLess< double >( 1.0 ) ) )
    {
      numDofs_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
        maxDofs_[ codim ] = 0;
      
      typedef typename LagrangePointSetContainerType::const_iterator IteratorType;
      const IteratorType end = lagrangePointSetContainer_.end();
      for( IteratorType it = lagrangePointSetContainer_.begin(); it != end; ++it )
      {
        const LagrangePointSetType *lagrangePointSet = *it;
        if( lagrangePointSet )
        {
          numDofs_ = std::max( numDofs_, lagrangePointSet->size() );
          for( int codim = 0; codim <= dimension; ++codim )
            maxDofs_[ codim ] = std::max( maxDofs_[ codim ], lagrangePointSet->maxDofs( codim ) );
        }
      }

      computeOffsets();
      dm_.addIndexSet( *this );
    }
    
    //! destructor 
    virtual ~LagrangeMapper ()
    {
      dm_.removeIndexSet( *this );
    }

    //! return overall number of degrees of freedom 
    int size () const
    {
      return size_;
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      // unsigned int codim, subEntity;
      const Fem::LocalKey &dofInfo = lagrangePointSet( entity.type() ).localKey( localDof );

      const unsigned int codim = dofInfo.codim();
      const int subIndex = indexSet_.subIndex( entity, dofInfo.subEntity(), codim );

      return offset_[ codim ] + subIndex * maxDofs_[ codim ] + dofInfo.index();
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const 
    {
      const unsigned int codim = Entity::codimension;
      assert( localDof < numEntityDofs( entity ) );
      return offset_[ codim ] + indexSet_.index( entity ) * maxDofs_[ codim ] + localDof;
    }
    
    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return numDofs_;
    }

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      return lagrangePointSet( entity.type() ).size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // This implementation only works for nonhybrid grids (simplices or cubes)
      return maxDofs_[ Entity::codimension ];
    }
    
    /** \brief Check, whether any DoFs are associated with a codimension */
    bool contains ( int codim ) const
    {
      return true;
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    bool fixedDataSize ( int codim ) const
    {
      return false;
    }

    /** \copydoc Dune::DofMapper::oldIndex */
    int oldIndex ( const int hole, const int codim ) const
    {
      const int n = maxDofs_[ codim ];
      return offset_[ codim ] + n*indexSet_.oldIndex( hole / n, codim ) + (hole % n);
    }

    /** \copydoc Dune::DofMapper::newIndex */
    int newIndex ( const int hole, const int codim ) const
    {
      const int n = maxDofs_[ codim ];
      return offset_[ codim ] + n*indexSet_.newIndex( hole / n, codim ) + (hole % n);
    }

    /** \copydoc Dune::DofMapper::numberOfHoles */
    int numberOfHoles ( const int codim ) const
    {
      return maxDofs_[ codim ] * indexSet_.numberOfHoles( codim );
    }

    /** \copydoc Dune::DofMapper::numBlocks
     */
    int numBlocks () const 
    {
      return dimension + 1;
    }

    /** \copydoc Dune::DofMapper::oldOffSet(const int block) const */
    int oldOffSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return oldOffSet_[ block ];
    }

    /** \copydoc Dune::DofMapper::offSet(const int block) const */
    int offSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return offset_[ block ];
    }

    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return BaseType::checkConsecutive( indexSet_ );
    }


    // Adaptation Methods (as for Index Sets)

    void insertEntity ( const EntityType &entity )
    {
      // check, whether we need to enlarge any block and compute oversized block size
      unsigned int oldUpperBound = size_;
      for( int codim = dimension ; codim >= 0; --codim )
      {
        const unsigned int codimSize = indexSet_.size( codim ) * maxDofs_[ codim ];
        const unsigned int newUpperBound = offset_[ codim ] + codimSize;
        if( newUpperBound > oldUpperBound )
          return computeOffsets( overShoot_ );
        oldUpperBound = offset_[ codim ];
      }
    }

    void removeEntity ( const EntityType &entity )
    {}

    void resize ()
    {
      computeOffsets();
    }

    bool compress ()
    {
      computeOffsets();
      return true;
    }

    void read_xdr ( const char *filename, int timestep )
    {
      computeOffsets();
    }

    void write_xdr ( const char *filename, int timestep )
    {}

  private:
    // prohibit copying and assignment
    LagrangeMapper ( const ThisType & );
    ThisType &operator=( const ThisType & );

    void computeOffsets ( const double overShoot = 1.0 )
    {
      // store old offset and calculate new offsets
      size_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        oldOffSet_[ codim ] = offset_[ codim ];
        offset_[ codim ] = size_;

        const unsigned int codimSize = indexSet_.size( codim ) * maxDofs_[ codim ];
        size_ += static_cast< unsigned int >( overShoot * codimSize );
      }
    }

    const LagrangePointSetType &lagrangePointSet ( const GeometryType type ) const
    {
      const LagrangePointSetType *lagrangePointSet
        = lagrangePointSetContainer_[ LocalGeometryTypeIndex::index( type ) ];
      assert( lagrangePointSet );
      return *lagrangePointSet;
    }

    // reference to dof manager needed for debug issues 
    DofManagerType &dm_;

    const IndexSetType &indexSet_;
    
    const LagrangePointSetContainerType &lagrangePointSetContainer_;

    // memory overshoot 
    const double overShoot_ ;

    unsigned int maxDofs_[ dimension+1 ];
    mutable unsigned int offset_[ dimension+1 ];
    mutable unsigned int oldOffSet_[ dimension+1 ];
    mutable unsigned int size_;
    unsigned int numDofs_;
  };
#endif // #ifdef USE_TWISTFREE_MAPPER

} // end namespace Dune 

#endif // #ifndef DUNE_LAGRANGESPACE_MAPPER_HH
