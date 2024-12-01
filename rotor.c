#include "rotor.h"
#include "mvec.h"
#include <math.h>

struct GaRotor newGaRotor(float angle, struct GaBivector rotPlane) {
  float rotNorm = Norm(rotPlane.mvec);

  struct GaRotor newRotor = {
      .mvec.mvec[0] = cosf(angle / 2.0),
      .mvec.mvec[4] = sinf(angle / 2.0) * (rotPlane.mvec.mvec[4] / rotNorm),
      .mvec.mvec[5] = sinf(angle / 2.0) * (rotPlane.mvec.mvec[5] / rotNorm),
      .mvec.mvec[6] = sinf(angle / 2.0) * (rotPlane.mvec.mvec[6] / rotNorm),
  };

  return newRotor;
}
