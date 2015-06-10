namespace Dune
{
  namespace Fem
  {

    template <class GridPart>
    typename CacheProvider<GridPart, 1>::MapperContainerType
    CacheProvider<GridPart, 1>::mappers_;


    template <class GridPart>
    typename CacheProvider<GridPart, 1>::MapperIteratorType
    CacheProvider<GridPart, 1>::createMapper(const QuadratureType& quad,
                                             GeometryType elementGeometry,
                                             std::integral_constant< bool, true > )
    {
      typedef TwistProvider<ct, dim-codim> TwistProviderType;
      typedef typename TwistProviderType::TwistStorageType TwistStorageType;

      const TwistStorageType& twistMappers =
        TwistProviderType::getTwistStorage(quad);
      const MapperVectorType pointMappers =
        PointProvider<ct, dim, codim>::getMappers(quad,
                                                  twistMappers.getPoints(),
                                                  elementGeometry);

      const int numFaces = pointMappers.size();
      const int maxTwist = twistMappers.maxTwist();
      const int minTwist = twistMappers.minTwist();

      QuadratureKeyType key ( elementGeometry, quad.id() );
      MapperIteratorType it = mappers_.insert
        (std::make_pair( key,
                         CacheStorageType(numFaces, maxTwist))).first;

      for (int face = 0; face < numFaces; ++face)
      {
        for (int twist = minTwist; twist < maxTwist; ++twist) {
          it->second.addMapper(pointMappers[face],
                               twistMappers.getMapper(twist),
                               face, twist);
        }
      }

      return it;
    }



    template <class GridPart>
    typename CacheProvider<GridPart, 1>::MapperIteratorType
    CacheProvider<GridPart, 1>::createMapper(const QuadratureType& quad,
                                             GeometryType elementGeometry,
                                             std::integral_constant< bool, false > )
    {
      const MapperVectorType pointMappers =
        PointProvider<ct, dim, codim>::getMappers(quad, elementGeometry);

      const int numFaces = pointMappers.size();

      QuadratureKeyType key ( elementGeometry, quad.id() );

      MapperIteratorType it
        = mappers_.insert(std::make_pair(key, CacheStorageType(numFaces))).first;

      for (int face = 0; face < numFaces; ++face)
        it->second.addMapper(pointMappers[face], face);

      return it;
    }

  } // namespace Fem

} // namespace Dune
