#ifndef DUNE_FEM_CACHINGQUADRATURE_HH
#define DUNE_FEM_CACHINGQUADRATURE_HH

//- Dune includes
#include <dune/common/misc.hh>

//- Local includes
#include "elementquadrature.hh"
#include "caching/twistutility.hh"
#include "caching/pointmapper.hh"
#include "caching/cacheprovider.hh"

#include "cachingpointlist.hh"

namespace Dune
{
  namespace Fem 
  {

    /** \class CachingQuadrature
     *  \ingroup Quadrature
     *  \brief quadrature class supporting base function caching
     *
     *  A CachingQuadrature is a conceptual extension to the ElementQuadrature.
     *  It provides an additional mapping from local quadrature point numbers on
     *  a subentity's reference element to global quadrature point numbers on the
     *  codim-0 reference element. Consider, for instance, a quadrature for one
     *  of the faces of a tetrahedron: It provides n local quadrature points, which
     *  can lie on one of the four faces, resulting in 4*n global quadrature points.
     *  
     *  The information from the mapping can be used to cache a base function on
     *  those global quadrature points.
     *
     *  \note If you don't want caching, you can use ElementQuadrature instead.
     *
     *  For the actual implementations, see
     *  - CachingQuadrature<GridPartImp,0>
     *  - CachingQuadrature<GridPartImp,1>
     */
    template< typename GridPartImp, int codim >
    class CachingQuadrature;



    /** \copydoc CachingQuadrature */
    template< typename GridPart >
    class CachingQuadrature< GridPart, 0 >
    : public CachingPointList< GridPart, 0, ElementQuadratureTraits< GridPart, 0 > >
    {
    public:
      //! type of grid partition
      typedef GridPart GridPartType;
      
      //! codimension of the element quadrature
      enum { codimension = 0 };
     
    private:
      typedef ElementQuadratureTraits< GridPartType, codimension > IntegrationTraits;
      
      typedef CachingQuadrature< GridPartType, codimension > ThisType;
      typedef CachingPointList< GridPartType, codimension, IntegrationTraits >
        BaseType;

    public:
      //! Dimension of the world.
      enum { dimension = BaseType::dimension };

      //! Just another name for double...
      typedef typename BaseType :: RealType RealType;
      //! The type of the coordinates in the codim-0 reference element.
      typedef typename BaseType :: CoordinateType CoordinateType;

      // for compatibility
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
      
    protected:
      using BaseType :: quadImp;

    public:
      /** \brief constructor
       *
       *  \param[in]  entity  entity, on whose reference element the quadratre
       *                      lives
       *  \param[in]  order   desired minimal order of the quadrature
       */
      CachingQuadrature( const EntityType &entity,
                         int order )
      : BaseType( entity.type(), order )
      {}

      /** \brief copy constructor
       *
       *  \param[in]  org  element quadrature to copy
       */
      CachingQuadrature( const ThisType &org ) 
      : BaseType( org )
      {}

      /** \copydoc Dune::Fem::ElementQuadrature<GridPartImp,0>::weight */
      const RealType &weight ( size_t i ) const 
      {
        return quadImp().weight( i );
      }
    };


   
    /** \copydoc CachingQuadrature */
    template< typename GridPartImp >
    class CachingQuadrature< GridPartImp, 1 >
    : public CachingPointList
      < GridPartImp, 1, ElementQuadratureTraits< GridPartImp, 1 > >
    {
    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the element quadrature
      enum { codimension = 1 };

    private:
      typedef ElementQuadratureTraits< GridPartType, codimension > IntegrationTraits;

      typedef CachingQuadrature< GridPartType, codimension > ThisType;
      typedef CachingPointList< GridPartType, codimension, IntegrationTraits >
        BaseType;

    protected:
      using BaseType :: quadImp;

    public:
      //! Dimeinsion of the world
      enum { dimension = BaseType::dimension };
      
      //! A double... or whatever your grid wants
      typedef typename BaseType::RealType RealType;

      //! The coordinates of the quadrature points in the codim-0 reference
      //! element
      typedef typename BaseType::CoordinateType CoordinateType;

      //! Type of the intersection iterator
      typedef typename BaseType :: IntersectionIteratorType IntersectionIteratorType;
      typedef typename IntersectionIteratorType :: Intersection IntersectionType;
      
      //! type of quadrature used for non-conforming intersections  
      typedef ElementQuadrature< GridPartImp, codimension > NonConformingQuadratureType; 

    public:
      /** \brief constructor
       *
       *  \note The CachingQuadrature requires the grid part to get twist
       *        information for TwistUtility (see also
       *        ElementQuadrature<GridPartImp,1>).
       * 
       *  \param[in]  gridPart      grid partition
       *  \param[in]  intersection  intersection
       *  \param[in]  order         desired order of the quadrature
       *  \param[in]  side          either INSIDE or OUTSIDE; codim-0 entity for 
       *                            which the ElementQuadrature shall be created
       */
      CachingQuadrature( const GridPartType &gridPart, 
                         const IntersectionType &intersection,
                         int order,
                         typename BaseType::Side side )
      : BaseType( gridPart, intersection, order, side )
      {}

      /** \brief copy constructor
       *
       *  \param[in]  org  element quadrature to copy
       */
      CachingQuadrature( const ThisType& org ) 
      : BaseType( org )
      {}

      /** \copydoc Dune::Fem::ElementQuadrature<GridPartImp,1>::weight */
      const RealType &weight( size_t i ) const 
      {
        return quadImp().weight(i);
      }
    };

  } //namespace Fem

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: CachingQuadrature ;
#endif // DUNE_FEM_COMPATIBILITY

} //namespace Dune

#endif // #ifndef DUNE_FEM_CACHINGQUADRATURE_HH
