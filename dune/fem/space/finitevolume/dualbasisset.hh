#ifndef DUNE_FEM_SPACE_FINITEVOLUME_DUALBASISSET_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_DUALBASISSET_HH

#include <cassert>
#include <cstddef>

#include <type_traits>
#include <utility>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/storage/entitygeometry.hh>
#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeBasisFunctionSet
    // ----------------------------

    template< class Entity, class Range >
    struct VertexCenteredFiniteVolumeBasisFunctionSet
      : public EntityGeometryStorage< Entity >
    {
    protected:
      typedef EntityGeometryStorage< Entity >  BaseType;

    public:
      /** \copydoc Dune::Fem::BasisFunctionSet::EntityType */
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::Geometry   Geometry;

      /** \copydoc Dune::Fem::BasisFunctionSet::FunctionSpaceType */
      typedef FunctionSpace< typename Entity::Geometry::ctype, typename Range::value_type,
                             Entity::Geometry::coorddimension, Range::dimension
                           > FunctionSpaceType;

      /** \copydoc Dune::Fem::BasisFunctionSet::DomainType */
      typedef typename FunctionSpaceType::DomainType DomainType;
      /** \copydoc Dune::Fem::BasisFunctionSet::RangeType */
      typedef typename FunctionSpaceType::RangeType RangeType;
      /** \copydoc Dune::Fem::BasisFunctionSet::JacobianRangeType */
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      /** \copydoc Dune::Fem::BasisFunctionSet::HessianRangeType */
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      /** \copydoc Dune::Fem::BasisFunctionSet::ReferenceElementType */
      typedef std::decay_t< decltype( Dune::ReferenceElements< typename EntityType::Geometry::ctype, EntityType::Geometry::coorddimension >::general( std::declval< const Dune::GeometryType & >() ) ) > ReferenceElementType;

    protected:
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
      typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

      int scalarSize_;
      int size_;

    public:
      /** \name Construction
       *  \{
       */

      VertexCenteredFiniteVolumeBasisFunctionSet () : scalarSize_(0), size_(0) {}

      explicit VertexCenteredFiniteVolumeBasisFunctionSet ( const EntityType &entity )
        : BaseType( entity ),
          scalarSize_( entity.subEntities( EntityType::dimension ) ),
          size_( scalarSize_ * RangeType::dimension )
      {}

      /** \} */

      /** \name Public member methods
       *  \{
       */
      using BaseType::entity;
      using BaseType::valid;
      using BaseType::geometry;
      using BaseType::referenceElement;
      using BaseType::type;

      /** \copydoc Dune::Fem::BasisFunctionSet::order */
      static constexpr int order () { return 0; }

      /** \copydoc Dune::Fem::BasisFunctionSet::size */
      std::size_t size () const { return size_; }

      /** \brief size of scalar basis function set */
      std::size_t scalarSize () const { return scalarSize_; }

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Quadrature, class Vector, class DofVector >
      void axpy ( const Quadrature &quadrature, const Vector &values, DofVector &dofs ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          axpyImpl( values[ qp ], dofs );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Quadrature, class VectorA, class VectorB, class DofVector >
      void axpy ( const Quadrature &quadrature, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          axpyImpl( valuesA[ qp ], dofs );
          axpyImpl( valuesB[ qp ], dofs );
        }
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        axpyImpl( valueFactor, dofs );
      }

    protected:
      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Point, class DofVector >
      void axpyImpl ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        // TODO: double check implementation
        assert( false );
        std::abort();

        for( int s = 0, i = 0; s < scalarSize_; ++s )
        {
          const RangeFieldType phi = shapeFunction( s, x );
          for( int r = 0; r < RangeType::dimension; ++r, ++i )
            dofs[ i ] += valueFactor[ r ] * phi ;
        }
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class DofVector >
      void axpyImpl ( const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {}

    public:
      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {}

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::evaluateAll */
      template< class Quadrature, class DofVector, class RangeArray >
      void evaluateAll ( const Quadrature &quadrature, const DofVector &dofs, RangeArray &ranges ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          evaluateAll( quadrature[ qp ], dofs, ranges[ qp ] );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::evaluateAll */
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        value = 0.;
        assert( dofs.size() == size_ );
        for( int i = 0, dof = 0; i < scalarSize_; ++i )
        {
          RangeFieldType phi = shapeFunction( i, x );
          for( int r=0; r<RangeType::dimension; ++r, ++dof  )
            value[ r ] += phi * dofs[ dof ];
        }
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::evaluateAll */
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        for( int s = 0, i = 0; s < scalarSize_; ++s )
        {
          RangeFieldType phi = shapeFunction( s, x );
          for( int r=0; r<RangeType::dimension; ++r, ++i )
          {
            values[ i ] = RangeType( 0 );
            values[ i ][ r ] = phi;
          }
        }
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::jacobianAll */
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quadrature, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          jacobianAll( quadrature[ qp ], dofs, jacobians[ qp ] );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::jacobianAll */
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( 0 );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::jacobianAll */
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          jacobians[ i ] = JacobianRangeType( 0 );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::hessianAll */
      template< class QuadratureType, class DofVector, class HessianArray >
      void hessianAll ( const QuadratureType &quadrature, const DofVector &dofs, HessianArray &hessians ) const
      {
        assert( hessians.size() >= quadrature.nop() );
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          hessians[qp] = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::hessianAll */
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        hessian = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::hessianAll */
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          hessians[ i ] = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }
    protected:
      template <class Point>
      RangeFieldType shapeFunction( const int phi, const Point& x ) const
      {
        if constexpr ( Entity::dimension == 2 )
        {
          if( entity().type().isSimplex() )
          {
            return triangleShapeFunction( phi, x );
          }
          else
          {
            return quadrilateralShapeFunction( phi, x );
          }
        }

        DUNE_THROW(NotImplemented,"shape functions not implemented yet!");
        return 0.0;
      }

      template <class Point, class AffineGeometry>
      RangeFieldType checkInsideTriangle(const Point& x, const AffineGeometry& tr_a, const AffineGeometry& tr_b ) const
      {
        const auto& refElem = referenceElement();
        auto xLoc = tr_a.local( Fem::coordinate( x ) );
        if( refElem.checkInside( xLoc ) )
        {
          //std::cout << "Inside tr_a " << std::endl;
          return 1.0;
        }

        xLoc = tr_b.local( Fem::coordinate( x ) );
        if( refElem.checkInside( xLoc ) )
        {
          //std::cout << "Inside tr_b " << std::endl;
          return 1.0;
        }

        return 0.0;
      }

      template <class Point>
      RangeFieldType triangleShapeFunction(const int phi, const Point& x ) const
      {
        static_assert( Entity::dimension == 2 );
        assert( entity().type().isSimplex() );

        typedef Dune::AffineGeometry< typename Geometry::ctype,
                Geometry::mydimension, Geometry::coorddimension > AffineGeometryType;

        typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

        static std::vector< AffineGeometryType > triangles;

        if( triangles.empty() )
        {
          GlobalCoordinate p0( 0 );
          GlobalCoordinate p1( {{1,0}} );
          GlobalCoordinate p2( {{0,1}} );
          GlobalCoordinate p3( {{0.5, 0}} );
          GlobalCoordinate p4( {{0.5, 0.5}} );
          GlobalCoordinate p5( {{0, 0.5}} );
          GlobalCoordinate p6( {{ 1./3. , 1./3.}} );
          ////////////////////////////////////////////////////////
          //  p_2
          //    *
          //    | *
          //    |   *
          //    |     *
          //    |T4     *
          //    |      T5 *
          // p_5+ .        .* p_4
          //    |   .p_6 .    *
          //    | T1   *    T3  *
          //    |        .        *
          //    |          .        *
          //    |      T0   . T2      *
          //    *-----------+-----------*
          //  p_0          p_3          p_1
          //
          /////////////////////////////////////////////////////////

          {
            std::vector< GlobalCoordinate > coords = {{ p0, p3, p6 }};
            triangles.emplace_back( AffineGeometryType( entity().type(), coords) );
          }
          {
            std::vector< GlobalCoordinate > coords = {{ p0, p6, p5 }};
            triangles.emplace_back( AffineGeometryType( entity().type(), coords) );
          }

          {
            std::vector< GlobalCoordinate > coords = {{  p1, p6, p3 }} ;
            triangles.emplace_back( AffineGeometryType( entity().type(), coords ));
          }
          {
            std::vector< GlobalCoordinate > coords = {{ p1, p4, p6 }};
            triangles.emplace_back( AffineGeometryType( entity().type(), coords ));
          }

          {
            std::vector< GlobalCoordinate > coords = {{ p2, p5, p6 }};
            triangles.emplace_back( AffineGeometryType( entity().type(), coords ));
          }
          {
            std::vector< GlobalCoordinate > coords = {{ p2, p6, p4 }};
            triangles.emplace_back( AffineGeometryType( entity().type(), coords ));
          }
        }

        assert( triangles.size() == 6 );
        return checkInsideTriangle(x, triangles[ phi*2 + 0 ], triangles[ phi*2 + 1 ] );
      }


      template <class Point >
      RangeFieldType checkInsideQuadrilateral(const Point& x, const GlobalCoordinate& min, const GlobalCoordinate& max ) const
      {
        const auto& refElem = referenceElement();
        const auto& point = Fem::coordinate( x );
        for( int i=0; i<GlobalCoordinate::dimension; ++i )
        {
          if( point[ i ] < min[ i ] || point[ i ] > max[ i ] )
            return 0.0;
        }

        return 1.0;
      }

      template <class Point>
      RangeFieldType quadrilateralShapeFunction(const int phi, const Point& x ) const
      {
        static_assert( Entity::dimension == 2 );
        assert( entity().type().isCube() );

        typedef Dune::AffineGeometry< typename Geometry::ctype,
                Geometry::mydimension, Geometry::coorddimension > AffineGeometryType;

        typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

        static std::vector< std::pair< GlobalCoordinate, GlobalCoordinate > > minMax;

        if( minMax.empty() )
        {
          ////////////////////////////////////////////////////////
          //  p_2                      p_3
          //    *-----------------------*
          //    |           |           |
          //    |           |           |
          //    |    Q2     |    Q3     |
          //    |           |           |
          //    |           |           |
          //    +-----------+-----------+
          //    |           |           |
          //    |           |           |
          //    |    Q0     |    Q1     |
          //    |           |           |
          //    |           |           |
          //    *-----------+-----------*
          //  p_0                       p_1
          //
          /////////////////////////////////////////////////////////

          // Q_0
          {
            GlobalCoordinate min( 0 );
            GlobalCoordinate max( 0.5 );
            minMax.push_back( std::make_pair( min, max ) );
          }

          // Q_1
          {
            GlobalCoordinate min( {{ 0.5, 0.0 }} );
            GlobalCoordinate max( {{ 1.0, 0.5 }} );
            minMax.push_back( std::make_pair( min, max ) );
          }

          // Q_2
          {
            GlobalCoordinate min( {{ 0.0, 0.5 }} );
            GlobalCoordinate max( {{ 0.5, 1.0 }} );
            minMax.push_back( std::make_pair( min, max ) );
          }

          // Q_3
          {
            GlobalCoordinate min( 0.5 );
            GlobalCoordinate max( 1.0 );
            minMax.push_back( std::make_pair( min, max ) );
          }
        }
        return checkInsideQuadrilateral(x, minMax[ phi ].first, minMax[ phi ].second );
      }

    };

    /** \} */

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSET_HH
