#include "vector.h"
#include "mvec.h"

struct GaVector newGaVec(float e1, float e2, float e3) {
  struct GaVector newVector = {
      .mvec.mvec[1] = e1,
      .mvec.mvec[2] = e2,
      .mvec.mvec[3] = e3,
  };

  return newVector;
}
