#ifndef VECTOR_H_
#define VECTOR_H_
#include "mvec.h"

struct GaVector {
  // There are 8 bases vectors in G3
  // 1, e1, e2, e3, e1e2, e3e1, e2e3, e1e2e3
  struct GaMultivector mvec;
};

struct GaVector newGaVec(float e1, float e2, float e3);

#endif // VECTOR_H_
