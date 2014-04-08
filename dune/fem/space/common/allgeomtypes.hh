#ifndef DUNE_FEM_ALLGEOMTYPES_HH
#define DUNE_FEM_ALLGEOMTYPES_HH

//- system includes 
#include <vector>
#include <map>

//- Dune includes 
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/capabilities.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int dim >
  class UGGrid;

  namespace Fem
  {

    // GeometryInformation
    // -------------------

    /**  \brief ReferenceVolume and local bary center keeper class. 
     */
    template< class GridImp, int codim >
    class GeometryInformation  
    {
      typedef GeometryInformation< GridImp, codim > ThisType;

    public:
      //! grid type 
      typedef GridImp GridType;

      //! dimension 
      static const int dim = GridType::dimension - codim;

      //! coordinate type 
      typedef typename GridType::ctype ctype;

      //! type of reference element 
      typedef Dune::ReferenceElement< ctype, dim > ReferenceElementType; 

      //! type of domain vector 
      typedef FieldVector<ctype, dim> DomainType;

    protected:
      //! constructor creating empty geometry information 
      GeometryInformation ()
      {}

    public:
      //! creating geometry information due to given geometry types list 
      explicit GeometryInformation( const std::vector< GeometryType > &geomTypes )
      {
        buildMaps( geomTypes );
      }

      //! return local bary center for geometry of type type 
      const DomainType &localCenter ( const GeometryType &type ) const 
      {
        return referenceElement( type ).position( 0, 0 );
      }

      //! return volume of reference element for geometry of type type 
      double referenceVolume ( const GeometryType &type ) const
      {
        return referenceElement( type ).volume();
      }

      //! return reference element for type 
      static const ReferenceElementType &referenceElement ( const GeometryType &type )
      {
        return Dune::ReferenceElements< ctype, dim >::general( type );
      }

    protected:  
      //! build maps 
      void buildMaps ( const std::vector< GeometryType > &geomTypes )
      {}
    };


    /**  \brief default implementation uses method geomTypes of given index
         set. Used in DiscreteFunctionSpaces.
     */
    template< class IndexSetImp, class GridImp >
    class AllGeomTypes : public GeometryInformation< GridImp , 0> 
    {
    public:
      typedef IndexSetImp IndexSetType;
      typedef GridImp GridType;

    private:
      typedef AllGeomTypes< IndexSetType, GridType > ThisType;

    protected:
      const IndexSetType &indexSet_;
      
    public:
      //! constructor storing index set reference 
      inline explicit AllGeomTypes( const IndexSetType &indexSet )
        : indexSet_( indexSet )
      {
        this->buildMaps( indexSet_.geomTypes(0) );
      }

      //! returns vector with geometry tpyes this index set has indices for
      const std :: vector< GeometryType > &geomTypes ( int codim ) const
      {
        return indexSet_.geomTypes( codim );
      }

      //! UGGrid might have different geom types 
      static bool multipleGeomTypes ()
      {
        assert( Dune::Capabilities :: hasSingleGeometryType < GridType > :: v );
        return false;
      }
    };


    /** \brief specialisation fir UGGrid, because geomTypes method of index
        sets not usable in this case. 
    */
    template< class IndexSetImp, int dim >
    class AllGeomTypes< IndexSetImp, Dune::UGGrid< dim > > 
    : public GeometryInformation< Dune::UGGrid< dim > , 0 > 
    {
      typedef AllGeomTypes< IndexSetImp, Dune::UGGrid< dim > > ThisType;

      static_assert( (dim == 2) || (dim == 3), "Invalid dimension for UG specified." );

    public:
      typedef IndexSetImp IndexSetType;
      typedef Dune::UGGrid< dim > GridType;
     
    protected:
      enum { ncodim = dim + 1 };
      
      std::vector< GeometryType > geomTypes_[ ncodim ];

    public:
      //! constructor storing index set reference 
      explicit AllGeomTypes ( const IndexSetType &indexSet )
      {
        // vertices 
        for (int i=0;i<dim;++i)
          geomTypes_[ i ].push_back( GeometryType( GeometryType::simplex, dim-i ) );

        if ( ! ( indexSet.geomTypes(0).size() == 1 && 
                 indexSet.geomTypes(0)[0].isSimplex() ) )
        {
          if( dim == 2 )
          {
            // elements 
            // geomTypes_[ 0 ].push_back( GeometryType( GeometryType::simplex, 2 ) );
            geomTypes_[ 0 ].push_back( GeometryType( GeometryType::cube, 2 ) );

            // faces 
            // geomTypes_[ 1 ].push_back( GeometryType( GeometryType::cube, 1 ) );
          }
          else if( dim == 3 )
          {
            // elements 
            // geomTypes_[ 0 ].push_back( GeometryType( GeometryType::simplex, 3 ) );
            geomTypes_[ 0 ].push_back( GeometryType( GeometryType::cube, 3 ) );
            geomTypes_[ 0 ].push_back( GeometryType( GeometryType::prism, 3 ) );
            geomTypes_[ 0 ].push_back( GeometryType( GeometryType::pyramid, 3 ) );

            // faces 
            // geomTypes_[ 1 ].push_back( GeometryType( GeometryType::simplex, 2 ) );
            geomTypes_[ 1 ].push_back( GeometryType( GeometryType::cube, 2 ) );

            // edges 
            // geomTypes_[ 2 ].push_back( GeometryType( GeometryType::cube, 1 ) );
          }
        }
        this->buildMaps( geomTypes_[ 0 ] );
      }                  
      
      //! returns vector with geometry types this index set has indices for
      const std::vector< GeometryType > &geomTypes ( int codim ) const
      {
        assert( (codim >= 0) && (codim < ncodim) );
        return geomTypes_[ codim ];
      }

      //! UGGrid might have different geom types 
      static bool multipleGeomTypes ()
      {
        return true;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ALLGEOMTYPES_HH
