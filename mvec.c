#include "mvec.h"
#include <math.h>

struct GaMultivector GeometricProduct(struct GaMultivector a,
                                      struct GaMultivector b) {
  struct GaMultivector res;
  res.mvec[0] = a.mvec[0] * b.mvec[0] // scalar/scalar
                + a.mvec[1] * b.mvec[1] + a.mvec[2] * b.mvec[2] +
                a.mvec[3] * b.mvec[3] // vector dot product
                - a.mvec[4] * b.mvec[4] - a.mvec[5] * b.mvec[5] -
                a.mvec[6] * b.mvec[6]    // bivector dot produt
                - a.mvec[7] * b.mvec[7]; // trivector/trivector
  res.mvec[1] = (a.mvec[0] * b.mvec[1] + a.mvec[1] * b.mvec[0])     // scalar/e1
                + (a.mvec[4] * b.mvec[2] - a.mvec[2] * b.mvec[4])   // e12e2
                + (a.mvec[3] * b.mvec[5] - a.mvec[5] * b.mvec[3])   // e3e13
                + (-a.mvec[6] * b.mvec[7] - a.mvec[7] * b.mvec[6]); // e23e123
  res.mvec[2] = (a.mvec[0] * b.mvec[2] + a.mvec[2] * b.mvec[0])     // scalar/e2
                + (a.mvec[1] * b.mvec[4] - a.mvec[4] * b.mvec[1])   // e1e12
                + (a.mvec[6] * b.mvec[3] - a.mvec[3] * b.mvec[6])   // e23e3
                + (-a.mvec[5] * b.mvec[7] - a.mvec[7] * b.mvec[5]); // e31e123
  res.mvec[3] = (a.mvec[0] * b.mvec[3] + a.mvec[3] * b.mvec[0])     // scalar/e3
                + (a.mvec[5] * b.mvec[1] - a.mvec[1] * b.mvec[5])   // e31e1
                + (a.mvec[2] * b.mvec[6] - a.mvec[6] * b.mvec[2])   // e2e23
                + (-a.mvec[4] * b.mvec[7] - a.mvec[7] * b.mvec[4]); // e12e123
  res.mvec[4] = (a.mvec[0] * b.mvec[4] + a.mvec[4] * b.mvec[0])    // scalar/e12
                + (a.mvec[1] * b.mvec[2] - a.mvec[2] * b.mvec[1])  // e1e2
                + (a.mvec[5] * b.mvec[6] - a.mvec[6] * b.mvec[5])  // e31e23
                + (a.mvec[3] * b.mvec[7] + a.mvec[7] * b.mvec[3]); // e3e123
  res.mvec[5] = (a.mvec[0] * b.mvec[5] + a.mvec[5] * b.mvec[0])    // scalar/e31
                + (a.mvec[3] * b.mvec[1] - a.mvec[1] * b.mvec[3])  // e3e1
                + (a.mvec[6] * b.mvec[4] - a.mvec[4] * b.mvec[6])  // e23e12
                + (a.mvec[2] * b.mvec[7] + a.mvec[7] * b.mvec[2]); // e2e123
  res.mvec[6] = (a.mvec[0] * b.mvec[6] + a.mvec[6] * b.mvec[0])    // scalar/e23
                + (a.mvec[2] * b.mvec[3] - a.mvec[3] * b.mvec[2])  // e2e3
                + (a.mvec[4] * b.mvec[5] - a.mvec[5] * b.mvec[4])  // e12e31
                + (a.mvec[1] * b.mvec[7] + a.mvec[7] * b.mvec[1]); // e1e123
  res.mvec[7] = a.mvec[7] * b.mvec[0] + a.mvec[6] * b.mvec[1] +
                a.mvec[5] * b.mvec[2] + a.mvec[4] * b.mvec[3] +
                a.mvec[3] * b.mvec[4] + a.mvec[2] * b.mvec[5] +
                a.mvec[1] * b.mvec[6] + a.mvec[0] * b.mvec[7];
  return res;
};

struct GaMultivector DotProduct(struct GaMultivector a,
                                struct GaMultivector b) {
  struct GaMultivector res;
  res.mvec[0] = a.mvec[0] * b.mvec[0] // scalar/scalar
                + a.mvec[1] * b.mvec[1] + a.mvec[2] * b.mvec[2] +
                a.mvec[3] * b.mvec[3] // vector dot product
                - a.mvec[4] * b.mvec[4] - a.mvec[5] * b.mvec[5] -
                a.mvec[6] * b.mvec[6]    // bivector dot produt
                - a.mvec[7] * b.mvec[7]; // trivector/trivector
  res.mvec[1] = (a.mvec[0] * b.mvec[1] + a.mvec[1] * b.mvec[0])     // scalar/e1
                + (a.mvec[4] * b.mvec[2] - a.mvec[2] * b.mvec[4])   // e12e2
                + (a.mvec[3] * b.mvec[5] - a.mvec[5] * b.mvec[3])   // e3e13
                + (-a.mvec[6] * b.mvec[7] - a.mvec[7] * b.mvec[6]); // e23e123
  res.mvec[2] = (a.mvec[0] * b.mvec[2] + a.mvec[2] * b.mvec[0])     // scalar/e2
                + (a.mvec[1] * b.mvec[4] - a.mvec[4] * b.mvec[1])   // e1e12
                + (a.mvec[6] * b.mvec[3] - a.mvec[3] * b.mvec[6])   // e23e3
                + (-a.mvec[5] * b.mvec[7] - a.mvec[7] * b.mvec[5]); // e31e123
  res.mvec[3] = (a.mvec[0] * b.mvec[3] + a.mvec[3] * b.mvec[0])     // scalar/e3
                + (a.mvec[5] * b.mvec[1] - a.mvec[1] * b.mvec[5])   // e31e1
                + (a.mvec[2] * b.mvec[6] - a.mvec[6] * b.mvec[2])   // e2e23
                + (-a.mvec[4] * b.mvec[7] - a.mvec[7] * b.mvec[4]); // e12e123
  res.mvec[4] =
      (a.mvec[0] * b.mvec[4] +
       a.mvec[4] *
           b.mvec[0]) // scalar/e12
                      // +(a.mvec[1]*b.mvec[2]-a.mvec[2]*b.mvec[1]) // e1e2
                      // +(a.mvec[5]*b.mvec[6]-a.mvec[6]*b.mvec[5]) // e31e23
      + (a.mvec[3] * b.mvec[7] + a.mvec[7] * b.mvec[3]); // e3e123
  res.mvec[5] =
      (a.mvec[0] * b.mvec[5] +
       a.mvec[5] *
           b.mvec[0]) // scalar/e31
                      // +(a.mvec[3]*b.mvec[1]-a.mvec[1]*b.mvec[3]) // e3e1
                      // +(a.mvec[6]*b.mvec[4]-a.mvec[4]*b.mvec[6]) // e23e12
      + (a.mvec[2] * b.mvec[7] + a.mvec[7] * b.mvec[2]); // e2e123
  res.mvec[6] =
      (a.mvec[0] * b.mvec[6] +
       a.mvec[6] *
           b.mvec[0]) // scalar/e23
                      // +(a.mvec[2]*b.mvec[3]-a.mvec[3]*b.mvec[2]) // e2e3
                      // +(a.mvec[4]*b.mvec[5]-a.mvec[5]*b.mvec[4]) // e12e31
      + (a.mvec[1] * b.mvec[7] + a.mvec[7] * b.mvec[1]); // e1e123
  res.mvec[7] = a.mvec[7] * b.mvec[0]
                // + a.mvec[6] * b.mvec[1]
                // + a.mvec[5] * b.mvec[2]
                // + a.mvec[4] * b.mvec[3]
                // + a.mvec[3] * b.mvec[4]
                // + a.mvec[2] * b.mvec[5]
                // + a.mvec[1] * b.mvec[6]
                + a.mvec[0] * b.mvec[7];

  return res;
};

struct GaMultivector Reverse(struct GaMultivector a) {
  struct GaMultivector res;
  res.mvec[0] = a.mvec[0];
  res.mvec[1] = a.mvec[1];
  res.mvec[2] = a.mvec[2];
  res.mvec[3] = a.mvec[3];
  res.mvec[4] = -a.mvec[4];
  res.mvec[5] = -a.mvec[5];
  res.mvec[6] = -a.mvec[6];
  res.mvec[7] = -a.mvec[7];
  return res;
};

// \[|A|^2=\left< A\^dag A \right>_0\]
float Norm(struct GaMultivector a) {
  float res = sqrtf(GeometricProduct(Reverse(a), a).mvec[0]);

  return res;
};

// \[A^{-1}=\frac{A^\dag}{|A|^2}\]
struct GaMultivector Inverse(struct GaMultivector a) {
  struct GaMultivector res;
  float squaredLength = 1.0 / GeometricProduct(Reverse(a), a).mvec[0];

  res.mvec[0] = Reverse(a).mvec[0] * squaredLength;
  res.mvec[1] = Reverse(a).mvec[1] * squaredLength;
  res.mvec[2] = Reverse(a).mvec[2] * squaredLength;
  res.mvec[3] = Reverse(a).mvec[3] * squaredLength;
  res.mvec[4] = Reverse(a).mvec[4] * squaredLength;
  res.mvec[5] = Reverse(a).mvec[5] * squaredLength;
  res.mvec[6] = Reverse(a).mvec[6] * squaredLength;
  res.mvec[7] = Reverse(a).mvec[7] * squaredLength;
  return res;
};
