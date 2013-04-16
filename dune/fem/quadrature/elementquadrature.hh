#ifndef DUNE_FEM_ELEMENTQUADRATURE_HH
#define DUNE_FEM_ELEMENTQUADRATURE_HH

#include "quadrature.hh"
#include "elementpointlist.hh"

namespace Dune
{

  namespace Fem
  { 

    /*! \class ElementQuadrature
     *  \ingroup Quadrature
     *  \brief quadrature on the codim-0 reference element
     *
     *  DUNE quadratures are defined per geometry type, using local coordinates
     *  for the quadrature points. To evaluate a base function in some quadrature
     *  point, the quadrature must return points within the codim-0 reference
     *  element.
     *
     *  Now, assume you want to integrate over the face of a tetrahedron. This
     *  means you need a quadrature for a triangle, but the quadrature points
     *  should be specified with respect to the tetrahedron, since we want to
     *  evaluate our function in these points. This is where the ElementQuadrature
     *  comes into play.
     *
     *  The ElementQuadrature takes a subentity and transforms the quadrature
     *  corresponding to the geometry to the codim-0 reference element.
     *  
     *  To achieve this goal, an element quadrature depends stronger on the
     *  context in which it is used. For example, for each face within a
     *  tetrahedron (though they are all the same) we need a different
     *  ElementQuadrature, since the coordinates of the quadrature points
     *  with respect to the codim-0 entity differ for each face.
     *
     *  \note Actually, codim-1 element quadratures depend on the intersection.
     *
     *  \note This quadrature does not support caching of base functions in
     *        quadrature points (see also CachingQuadrature).
     *
     *  For the actual implementations, see
     *  - ElementQuadrature<GridPartImp,0>
     *  - ElementQuadrature<GridPartImp,1>
     */
    template< typename GridPartImp, int codim >
    class ElementQuadrature;



    template< class GridPartImp, int codim >
    struct ElementQuadratureTraits
    {
      // type of single coordinate
      typedef typename GridPartImp :: ctype ctype;

      // dimension of quadrature 
      enum { dimension = GridPartImp ::  dimension };

      // codimension of quadrature
      enum { codimension = codim };

      // type of used integration point list
      typedef Quadrature< ctype, dimension-codim, DefaultQuadratureTraits > IntegrationPointListType;
      
      // type of local coordinate (with respect to the codim-0 entity)
      typedef typename Quadrature< ctype, dimension, DefaultQuadratureTraits > :: CoordinateType
        CoordinateType; 
    };

    

    /** \copydoc ElementQuadrature */
    template< typename GridPartImp >
    class ElementQuadrature< GridPartImp, 0 >
    : public ElementIntegrationPointList< GridPartImp, 0, ElementQuadratureTraits< GridPartImp, 0 > >
    {
      typedef ElementQuadrature< GridPartImp, 0 > ThisType;
      typedef ElementQuadratureTraits< GridPartImp, 0 > IntegrationTraits;
      typedef ElementIntegrationPointList< GridPartImp, 0, IntegrationTraits > BaseType;

    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the element quadrature
      enum { codimension = 0 };

      //! dimension of the world
      enum { dimension = GridPartType :: dimension };
      
      //! type for reals (usually double)
      typedef typename GridPartType :: ctype RealType;
      
      //! type for coordinates in the codim-0 reference element 
      typedef typename IntegrationTraits :: CoordinateType CoordinateType;
      
      // for compatibility
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

    protected:
      using BaseType::quadImp;

    public:
      /*! \brief constructor
       *
       *  \param[in]  entity  entity, on whose reference element the quadratre
       *                      lives
       *  \param[in]  order   desired minimal order of the quadrature
       */
      ElementQuadrature( const EntityType &entity, int order )
      : BaseType( entity.type(), order )
      {}
      
      /** \brief copy constructor
       *
       *  \param[in]  org  element quadrature to copy
       */
      ElementQuadrature( const ThisType &org )
      : BaseType( org )
      {}

      /** \copydoc Dune::Fem::Quadrature::weight */
      const RealType &weight( size_t i ) const
      {
        return quadImp().weight( i );
      }
    };



    /** \copydoc ElementQuadrature */
    template< class GridPartImp >
    class ElementQuadrature< GridPartImp, 1 >
    : public ElementIntegrationPointList
      < GridPartImp, 1, ElementQuadratureTraits< GridPartImp, 1 > >
    {
    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the quadrature
      enum { codimension = 1 };

    private:
      typedef ElementQuadratureTraits< GridPartType, codimension > IntegrationTraits;
      
      typedef ElementQuadrature< GridPartType, codimension > ThisType;
      typedef ElementIntegrationPointList< GridPartType, 1, IntegrationTraits >
        BaseType;

    protected:
      using BaseType :: quadImp;
      
    public:
      //! dimension of the world
      enum { dimension = GridPartType :: dimension };
      
      //! type for reals (usually double)
      typedef typename GridPartType :: ctype RealType;

      //! type of the intersection iterator
      typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;

      //! type of coordinates in codim-0 reference element
      typedef typename IntegrationTraits :: CoordinateType CoordinateType;

      //! type of coordinate in codim-1 reference element
      typedef typename IntegrationTraits :: IntegrationPointListType :: CoordinateType
        LocalCoordinateType;

      //! type of quadrature for use on non-conforming intersections 
      typedef ThisType NonConformingQuadratureType;
      
    public:
      /*! \brief constructor
       *
       *  \param[in]  gridPart      grid partition (a dummy here)
       *  \param[in]  intersection  intersection
       *  \param[in]  order         desired order of the quadrature
       *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for 
       *                            which the ElementQuadrature shall be created
       */
      ElementQuadrature ( const GridPartType &gridPart,
                          const IntersectionType &intersection, 
                          int order, 
                          typename BaseType :: Side side )
      : BaseType( gridPart, intersection, order, side )
      {}

      /*! \brief copy constructor
       *
       *  \param[in]  org  element quadrature to copy
       */
      ElementQuadrature( const ElementQuadrature &org )
      : BaseType( org )
      {
      }
      
      /*! obtain the weight of the i-th quadrature point
       *
       *  \note The quadrature weights sum up to the volume of the corresponding
       *        reference element.
       *
       *  \param[in]  i  index of the quadrature point
       *
       *  \returns weight of the i-th quadrature point within the quadrature
       */
      const RealType &weight( size_t i ) const
      {
        return quadImp().weight( i );
      }
    };

  } // namespace Fem

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 
  
  using Fem :: ElementQuadrature ;
#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune

#endif // #ifndef DUNE_FEM_ELEMENTQUADRATURE_HH
