#ifndef ROTOR_H_
#define ROTOR_H_
#include "bivector.h"
#include "mvec.h"

struct GaRotor {
  // There are 8 bases vectors in G3
  // 1, e1, e2, e3, e1e2, e3e1, e2e3, e1e2e3
  struct GaMultivector mvec;
};

struct GaRotor newGaRotor(float angle, struct GaBivector rotPlane);
#endif // ROTOR_H_
