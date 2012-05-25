#ifndef DUNE_LAGRANGESPACE_RESTRICTPROLONG_HH
#define DUNE_LAGRANGESPACE_RESTRICTPROLONG_HH

#include <map>

#if HAVE_DUNE_GEOMETRY
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#else
#include <dune/common/geometrytype.hh>
#include <dune/grid/common/genericreferenceelements.hh>
#endif

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/function/localfunction/localfunction.hh>
#include <dune/fem/space/lagrangespace/lagrangepoints.hh>

namespace Dune
{

  template< class G, int ord >
  struct LagrangeLocalRestrictProlong
  {
    typedef G GridType;

    typedef typename GridType::ctype ctype;
    static const int dimension = GridType::dimension;

    typedef FieldVector< ctype, dimension > DomainVector;

    typedef Dune::LagrangePointSet< LeafGridPart< GridType >, ord > LagrangePointSet;

  private:
    typedef typename LagrangePointSet::template Codim< 0 >::SubEntityIteratorType
      EntityDofIterator;

    typedef std::map< const GeometryType, const LagrangePointSet * > LagrangePointSetMap;

  public:
    ~LagrangeLocalRestrictProlong ()
    {
      typedef typename LagrangePointSetMap::iterator Iterator;
      const Iterator end = lagrangePointSet_.end();
      for( Iterator it = lagrangePointSet_.begin(); it != end; ++it )
        delete it->second;
    }

    template< class DomainField >
    void setFatherChildWeight ( const DomainField &weight ) {}

    template< class FT, class ST, class LocalGeometry >
    void restrictLocal ( LocalFunction< FT > &lfFather,
                         const LocalFunction< ST > &lfSon,
                         const LocalGeometry &geometryInFather,
                         const bool initialize ) const
    {
      static const int dimRange = LocalFunction< ST >::dimRange;

      const GenericReferenceElement< ctype, dimension > &refSon
        = GenericReferenceElements< ctype, dimension >::general( lfSon.entity().type() );

      const LagrangePointSet &pointSet = lagrangePointSet( lfFather.entity() );

      const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
      for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
      {
        const unsigned int dof = *sit;
        const DomainVector &pointInFather = pointSet.point( dof );
        const DomainVector pointInSon = geometryInFather.local( pointInFather );
        if( refSon.checkInside( pointInSon ) )
        {
          typename LocalFunction< ST >::RangeType phi;
          lfSon.evaluate( pointInSon, phi );
          for( int coordinate = 0; coordinate < dimRange; ++coordinate )
            lfFather[ dimRange * dof + coordinate ] = phi[ coordinate ];
        }
      }
    }

    template< class FT, class ST, class LocalGeometry >
    void prolongLocal ( const LocalFunction< FT > &lfFather, LocalFunction< ST > &lfSon,
                        const LocalGeometry &geometryInFather,
                        bool initialize ) const
    {
      static const int dimRange = LocalFunction< FT >::dimRange;

      const LagrangePointSet &pointSet = lagrangePointSet( lfSon.entity() );

      const EntityDofIterator send = pointSet.template endSubEntity< 0 >( 0 );
      for( EntityDofIterator sit = pointSet.template beginSubEntity< 0 >( 0 ); sit != send; ++sit )
      {
        const unsigned int dof = *sit;
        const DomainVector &pointInSon = pointSet.point( dof );
        const DomainVector pointInFather = geometryInFather.global( pointInSon );
        
        typename LocalFunction< FT >::RangeType phi;
        lfFather.evaluate( pointInFather, phi );
        for( int coordinate = 0; coordinate < dimRange; ++coordinate )
          lfSon[ dimRange * dof + coordinate ] = phi[ coordinate ];
      }
    }

    bool needCommunication () const { return false; }

  protected:
    template< class Entity >
    const LagrangePointSet &lagrangePointSet ( const Entity &entity ) const
    {
      return lagrangePointSet( entity.type() );
    }

    const LagrangePointSet &lagrangePointSet ( const GeometryType &type ) const
    {
      typedef typename LagrangePointSetMap::iterator Iterator;
      Iterator it = lagrangePointSet_.find( type );
      if( it == lagrangePointSet_.end() )
        it = lagrangePointSet_.insert( it, std::make_pair( type, new LagrangePointSet( type, ord ) ) );
      assert( it->second != 0 );
      return *(it->second);
    }

  private:
    mutable LagrangePointSetMap lagrangePointSet_;
  };

} // namespace Dune

#endif // #ifndef DUNE_LAGRANGESPACE_RESTRICTPROLONG_HH
