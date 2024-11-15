#ifndef BIVECTOR_H_
#define BIVECTOR_H_
#include "mvec.h"

struct GaBivector {
  // There are 8 bases vectors in G3
  // 1, e1, e2, e3, e1e2, e3e1, e2e3, e1e2e3
  struct GaMultivector mvec;
};

struct GaBivector newGaBivec(float e1e2, float e3e1, float e2e3);

#endif // BIVECTOR_H_
