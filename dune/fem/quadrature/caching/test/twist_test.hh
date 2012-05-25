#ifndef DUNE_TWIST_TEST_HH
#define DUNE_TWIST_TEST_HH

#include <dune/fem/misc/test.hh>
#include "../twistprovider.hh"

namespace Dune {
  namespace Fem {

  class TwistProvider_Test : public Test
  {
  public:
    virtual void run();
   
    void lineTest();
    void triangleTest();
    void quadrilateralTest();
    void nonSymmetricTest();

  private:
  };

  } // end namespace Fem
} // end namespace Dune

#endif
