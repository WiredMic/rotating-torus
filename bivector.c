#include "bivector.h"
#include "mvec.h"

struct GaBivector newGaBivec(float e1e2, float e3e1, float e2e3) {
  struct GaBivector newBivector = {
      .mvec.mvec[4] = e1e2,
      .mvec.mvec[5] = e3e1,
      .mvec.mvec[6] = e2e3,
  };

  return newBivector;
}
