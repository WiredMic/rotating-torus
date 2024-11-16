#ifndef BIVECTOR_H_
#define BIVECTOR_H_
#include "mvec.h"
#include "vector.h"

struct GaBivector {
  // There are 8 bases vectors in G3
  // 1, e1, e2, e3, e1e2, e3e1, e2e3, e1e2e3
  struct GaMultivector mvec;
};

struct GaBivector newGaBivec(float e1e2, float e3e1, float e2e3);

struct GaBivector newGaBivec2(struct GaVector vec1, struct GaVector vec2);

#endif // BIVECTOR_H_
