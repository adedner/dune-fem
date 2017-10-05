#ifndef DUNE_FEM_IDPROVIDER_HH
#define DUNE_FEM_IDPROVIDER_HH

#include <cstdlib>

#include <dune/common/visibility.hh>

namespace Dune
{

  namespace Fem
  {

    //! Singleton that manages a globally unique identifier.
    class IdProvider
    {
    public:
      //! Access to the singleton object.
      DUNE_EXPORT static IdProvider& instance()
      {
        static IdProvider idProvider;
        return idProvider;
      }

      //! Return a new identifier.
      //! \note Identifiers are never freed.
      size_t newId() { return lowestFreeId_++; }

    private:
      //! Constructor (for the singleton object)
      IdProvider() :
        lowestFreeId_(0)
      {}

      IdProvider(const IdProvider&);
      IdProvider& operator=(const IdProvider&);

    private:
      size_t lowestFreeId_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IDPROVIDER_HH
