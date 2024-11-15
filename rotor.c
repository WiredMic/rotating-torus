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

  /* j: 5 */
  /* index: 137 */
  /*   e1: 2.754566, e2: 0.148581, e3: 0.470182  */
  /* j: 6 */
  /* index: 138 */
  /*   e1: 2.686310, e2: 0.000000, e3: 0.458531  */

  return newRotor;

  struct GaRotor newRotor1;
  newRotor.mvec.mvec[0] = cosf(angle / 2.0);
  newRotor.mvec.mvec[4] = sinf(angle / 2.0) * (rotPlane.mvec.mvec[4] / rotNorm);
  newRotor.mvec.mvec[5] = sinf(angle / 2.0) * (rotPlane.mvec.mvec[5] / rotNorm);
  newRotor.mvec.mvec[6] = sinf(angle / 2.0) * (rotPlane.mvec.mvec[6] / rotNorm);

  // get wrong output
  /* index: 138 */
  /*   e1: 260.894653, e2: -0.000001, e3: 2.657002  */
  /* j: 7 */
  /* index: 139 */
  /*   e1: 279.397614, e2: 2.159741, e3: 3.914960  */
  /* j: 8 */

  // or
  /* j: 6 */
  /* index: 138 */
  /*   e1: -nan, e2: nan, e3: -nan  */
  /* j: 7 */
  /* index: 139 */
  /*   e1: -nan, e2: nan, e3: -nan  */
  /* j: 8 */

  return newRotor1;
}
