#include <dune/fem/quadrature/caching/pointprovider.hh>
#include <dune/fem/misc/threadmanager.hh>

namespace Dune
{

  namespace Fem
  {

    //- class SimpleStorage
    template <class FunctionSpaceImp>
    StorageBase<FunctionSpaceImp>::StorageBase(const FactoryType& factory) :
      storage_( factory.numBaseFunctions() ),
      elementGeometry_(factory.geometry()),
      rangeTmp_( Fem :: ThreadManager :: maxThreads() ),
      jacobianTmp_( Fem :: ThreadManager :: maxThreads() )
    {
      for (int i = 0; i < factory.numBaseFunctions(); ++i) 
      {
        storage_[ i ] = factory.baseFunction(i);
      }

      // initialize quad ids with non-valid numbers
      for( size_t i = 0; i<rangeTmp_.size(); ++i )
      {
        rangeTmp_[ i ].second = ~0u ;
        jacobianTmp_[ i ].second = ~0u ;
      }
    }

    template <class FunctionSpaceImp>
    StorageBase<FunctionSpaceImp>::~StorageBase()
    {
      for (size_t i = 0; i <storage_.size(); ++i) 
      {
        // delete base functions 
        delete storage_[ i ];
        storage_[ i ] = 0;
      }
    }

    template <class FunctionSpaceImp>
    int StorageBase<FunctionSpaceImp>::numBaseFunctions() const 
    {
      return storage_.size();
    }

    template <class FunctionSpaceImp>
    template <int diffOrd>
    void StorageBase<FunctionSpaceImp>::
    evaluate(int baseFunct, const FieldVector<int, diffOrd>& diffVar, 
             const DomainType& xLocal, RangeType& result) const
    {
      assert(baseFunct >= 0 && baseFunct < numBaseFunctions());
      storage_[baseFunct]->evaluate(diffVar, xLocal, result);
    }

    template <class FunctionSpaceImp>
    void StorageBase<FunctionSpaceImp>::
    jacobian(int baseFunct, const DomainType& xLocal, 
             JacobianRangeType& result) const
    {
      assert(baseFunct >= 0 && baseFunct < numBaseFunctions());
      RangeType tmp;
      
      FieldVector<int, 1> diffVar1( 0 );
      for (int i = 0; i < DomainType::dimension; ++i) {
        diffVar1[0] = i;
        storage_[baseFunct]->evaluate(diffVar1, xLocal, tmp);
        for (int j = 0; j < RangeType::dimension; ++j) {
          result[j][i] = tmp[j];
        }
      }
    }


    
    // caching storage
    // ---------------

    template< class FunctionSpaceImp >
    template< class QuadratureType >
    inline void CachingStorage< FunctionSpaceImp >
      :: evaluate( const int baseFunction,
                   const FieldVector< int, 0 > &diffVariable,
                   const QuadraturePointWrapper< QuadratureType > &x,
                   RangeType &ret ) const
    {
      enum { cachable = Conversion< QuadratureType, Fem::CachingInterface > :: exists };
      const QuadratureType &quad = x.quadrature();
      const int pt = x.point();

      assert( !cachable || (rangestored_.find( quad.id() ) != rangestored_.end()) );
      Evaluate< QuadratureType, cachable >
        :: evaluate( *this, baseFunction, diffVariable, quad, pt, ranges_, ret );
    }
    

    
    template< class FunctionSpaceImp >
    template< class QuadratureType >
    inline void CachingStorage< FunctionSpaceImp >
      :: evaluate( const int baseFunction,
                   const FieldVector< int, 1 > &diffVariable,
                   const QuadraturePointWrapper< QuadratureType > &x,
                   RangeType &ret ) const
    {
      enum { cachable = Conversion< QuadratureType, Fem::CachingInterface > :: exists };
      const QuadratureType &quad = x.quadrature();
      const int pt = x.point();

      assert( !cachable || (jacobianstored_.find( quad.id() ) != jacobianstored_.end()) );
      Evaluate< QuadratureType, cachable >
        :: evaluate( *this, baseFunction, diffVariable, quad, pt, jacobians_, ret );
    }


    template< class FunctionSpaceImp >
    template< int diffOrder, class QuadratureType >
    inline void CachingStorage< FunctionSpaceImp >
      :: evaluate( const int baseFunction,
                   const FieldVector< int, diffOrder > &diffVariable,
                   const QuadraturePointWrapper< QuadratureType > &x,
                   RangeType &ret ) const
    {
      evaluate( baseFunction, diffVariable, coordinate( x ), ret );
    }


    template< class FunctionSpaceImp >
    template< class QuadratureType >
    inline void CachingStorage< FunctionSpaceImp >
      :: jacobian( const int baseFunction,
                   const QuadraturePointWrapper< QuadratureType > &x,
                   JacobianRangeType &ret ) const
    {
      enum { cachable = Conversion< QuadratureType, Fem::CachingInterface > :: exists };
      const QuadratureType &quad = x.quadrature();
      const int pt = x.point();

      assert( !cachable || (jacobianstored_.find( quad.id() ) != jacobianstored_.end()) );
      Evaluate<QuadratureType, cachable >
        :: jacobian( *this, baseFunction, quad, pt, jacobians_, ret );
    }

    

    template <class FunctionSpaceImp>
    inline void CachingStorage<FunctionSpaceImp>
      ::cacheQuadrature ( std::size_t id, std::size_t codim, std::size_t quadSize )
    {
      RangeIteratorType it = rangestored_.find(id);
      if (it == rangestored_.end()) 
      {
        it = addEntryInterface(id, codim, quadSize).first;
      }

      assert(rangestored_.find(id) != rangestored_.end());
      assert(jacobianstored_.find(id) != jacobianstored_.end());
    }

    template <class FunctionSpaceImp>
    inline typename CachingStorage<FunctionSpaceImp>::ReturnPairType
    CachingStorage<FunctionSpaceImp>::
    addEntryInterface(const size_t id, const size_t codim, const size_t quadSize ) const 
    {
      enum { dimension = DomainType::dimension };
      switch (codim) 
      {
        case 0: return addEntry<0>(id, quadSize);
        case 1: return addEntry<1>(id, quadSize);
        default: assert(false); abort();
      }

      assert(false);
      abort();
      // only fake 
      return addEntry<0>(id, quadSize);
    }

    //--addEntry
    template <class FunctionSpaceImp>
    template <int codimension>
    inline typename CachingStorage<FunctionSpaceImp>::ReturnPairType
    CachingStorage<FunctionSpaceImp>::
    addEntry(const size_t quadId, const size_t quadSize) const 
    {
      enum { dimension = DomainType::dimension };

      typedef typename FunctionSpaceImp :: DomainFieldType DomainFieldType;
      typedef Fem::PointProvider<DomainFieldType, dimension, codimension> PointProviderType;
      typedef typename PointProviderType::GlobalPointVectorType PointVectorType;

      const PointVectorType& points = 
        PointProviderType::getPoints(quadId, this->elementGeometry_);

      assert(rangestored_.find(quadId) == rangestored_.end());
      RangeIteratorType rit =
        rangestored_.insert(std::make_pair(quadId,true)).first; 
      assert(rangestored_.find(quadId) != rangestored_.end());

      assert(jacobianstored_.find(quadId) == jacobianstored_.end());
      JacobianRangeIteratorType jit =
        jacobianstored_.insert(std::make_pair(quadId,true)).first;

      assert(jacobianstored_.find(quadId) != jacobianstored_.end());

      FieldVector<int, 0> diffVar;

      if ((size_t)ranges_.size()   <= quadId) ranges_.resize(quadId+10);
      if ((size_t)jacobians_.size()<= quadId) jacobians_.resize(quadId+10);

      ranges_[quadId].resize(points.size());
      jacobians_[quadId].resize(points.size());

      const size_t pointSize = points.size();
      const size_t numBaseFunctions = this->numBaseFunctions(); 
      for (size_t i = 0; i < pointSize; ++i) 
      {
        ranges_[quadId][i].resize(numBaseFunctions);
        jacobians_[quadId][i].resize(numBaseFunctions);

        for (size_t j = 0; j < numBaseFunctions; ++j) 
        {
          // evaluate value and jacobian and store it 
          this->evaluate(j, diffVar, points[i], ranges_[quadId][i][j]);
          this->jacobian(j, points[i], jacobians_[quadId][i][j]);
        }
      }

      return std::make_pair(rit, jit);
    }

  } // namespace Fem
  
} // namespace Dune