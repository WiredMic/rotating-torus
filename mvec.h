#ifndef MVEC_H_
#define MVEC_H_

struct GaMultivector {
  // There are 8 bases vectors in G3
  // 1, e1, e2, e3, e1e2, e3e1, e2e3, e1e2e3
  float mvec[8];
};

struct GaMultivector GeometricProduct(struct GaMultivector a,
                                      struct GaMultivector b);

struct GaMultivector DotProduct(struct GaMultivector a,
                                struct GaMultivector b);
struct GaMultivector Reverse(struct GaMultivector a);
float Norm(struct GaMultivector a);
struct GaMultivector Inverse(struct GaMultivector a);
#endif // MVEC_H_
