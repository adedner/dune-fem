#ifndef DUNE_POINTMAPPER_HH
#define DUNE_POINTMAPPER_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/fem/quadrature/quadrature.hh>

namespace Dune {

  class QuadratureKey  
  {
  protected:
    const size_t id_;
    
    static inline size_t quadId2MapperId(const GeometryType& elemGeo, const size_t quadId) 
    {
      // split between different geometry basic types 
      assert( (int)elemGeo.basicType() >= 0 );
      assert( quadId < 65536 );
      return elemGeo.basicType() * 65536 + quadId; 
    }
  public:  
  
    QuadratureKey(const GeometryType& geoType, const size_t id) 
      : id_( quadId2MapperId(geoType, id)) 
    {}
    
    QuadratureKey(const QuadratureKey& other) 
      : id_(other.id_)
    {}
    
    bool operator < (const QuadratureKey& other) const 
    {
      return (id_ < other.id_);
    }
    
    bool operator == (const QuadratureKey& other) const 
    {
      return id_ == other.id_;
    }
  };
    
  template <class ct, int dim>
  struct CachingTraits {
    //! type of integration point list implementation, fix type here 
    typedef IntegrationPointListImp<ct, dim> QuadratureType;
    //! extracted types from integration point list 
    typedef typename QuadratureType::CoordinateType PointType;
    typedef std::vector<PointType>    PointVectorType;
    typedef std::vector<size_t>       MapperType;
    typedef std::vector<MapperType>   MapperVectorType;

    typedef QuadratureKey QuadratureKeyType;

    // minimal twist is -4 for hexahedrons 
    // so we add 4 to start from zero 
    enum { twistOffset_ = 4 };

  };
  

} // end namespace Dune

#endif
