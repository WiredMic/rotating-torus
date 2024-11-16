#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

// Include basic geometric algebra structs
#include "bivector.h"
#include "mvec.h"
#include "rotor.h"
#include "vector.h"

// define constants
#define TAU 6.283185
#define minorRadius 1.0
#define majorRadius 3.0
#define NUM_POINT_IN_CIRC 12
#define NUM_CIRC_IN_TORUS 12
#define NUM_POINT_IN_TORUS (NUM_POINT_IN_CIRC * NUM_CIRC_IN_TORUS)
#define DONUT_WIDTH 25
#define DONUT_WIDTH_BUFFER 10
#define DONUT_HEIGTH (DONUT_WIDTH / 2.0)

// pass by reference
// \[\vec{v}' = R_2^\dag(R_1^\dag \vec{v} R_1 + \vec{r})R_2\]
// \[ = (R_1R_2)^\dag \vec{v} R_1R_2 + R_2^\dag\vec{r}R_2\]
void makeTorus(struct GaVector *pointArray) {

  struct GaVector firstPoint = newGaVec(minorRadius, 0.0, 0.0);
  struct GaVector translate = newGaVec(majorRadius, 0.0, 0.0);

  float angle1 = TAU / NUM_POINT_IN_CIRC;
  struct GaBivector rotPlane1 = newGaBivec(1.0, 0.0, 0.0);
  /* printf("rotPlane1: e12: %f, e31: %f, e23: %f \n", rotPlane1.mvec.mvec[4],
   */
  /*        rotPlane1.mvec.mvec[5], rotPlane1.mvec.mvec[6]); */
  float angle2 = TAU / NUM_CIRC_IN_TORUS;
  struct GaBivector rotPlane2 = newGaBivec(0.0, 1.0, 0.0);
  /* printf("rotPlane2: e12: %f, e31: %f, e23: %f \n", rotPlane2.mvec.mvec[4],
   */
  /*        rotPlane2.mvec.mvec[5], rotPlane2.mvec.mvec[6]); */

  for (int i = 0; i < NUM_CIRC_IN_TORUS; i++) {
    // rotation around in the e3e1 plane
    struct GaRotor rotor2 = newGaRotor(angle2 * i, rotPlane2);

    /* printf("    s: %f, e1: %f, e2: %f, e3: %f, e12: %f, e31: %f, e23: %f \n",
     */
    /*        rotor2.mvec.mvec[0], rotor2.mvec.mvec[1], rotor2.mvec.mvec[2], */
    /*        rotor2.mvec.mvec[3], rotor2.mvec.mvec[4], rotor2.mvec.mvec[5], */
    /*        rotor2.mvec.mvec[6]); */

    // R_2^\dag\vec{r}R_2
    struct GaVector translateRot = {
        .mvec = GeometricProduct(
            GeometricProduct(Reverse(rotor2.mvec), translate.mvec),
            rotor2.mvec),
    };

    /* printf("i: %d\n", i); */

    for (int j = 0; j < NUM_POINT_IN_CIRC; j++) {
      // rotation around the
      struct GaRotor rotor1 = newGaRotor(angle1 * j, rotPlane1);

      // The rotor that rotates from first point to the point in the torus
      // \[R_1R_2\]
      struct GaRotor rotor12 = {
          .mvec = GeometricProduct(rotor1.mvec, rotor2.mvec),
      };

      // \[(R_1R_2)^\dag \vec{v} R_1R_2\]
      struct GaVector vectorRot = {
          .mvec = GeometricProduct(
              GeometricProduct(Reverse(rotor12.mvec), firstPoint.mvec),
              rotor12.mvec),
      };

      // \[ = (R_1R_2)^\dag \vec{v} R_1R_2 + R_2^\dag\vec{r}R_2\]
      struct GaVector pointVec = {
          .mvec.mvec[1] = vectorRot.mvec.mvec[1] + translateRot.mvec.mvec[1],
          .mvec.mvec[2] = vectorRot.mvec.mvec[2] + translateRot.mvec.mvec[2],
          .mvec.mvec[3] = vectorRot.mvec.mvec[3] + translateRot.mvec.mvec[3],
      };

      pointArray[j + i * NUM_POINT_IN_CIRC] = pointVec;

      /* printf("  j: %d\n", j); */
      /* printf("  index: %d\n", j + i * NUM_POINT_IN_CIRC); */
      /* printf("    e1: %f, e2: %f, e3: %f \n", */
      /*        pointArray[j + i * NUM_POINT_IN_CIRC].mvec.mvec[1], */
      /*        pointArray[j + i * NUM_POINT_IN_CIRC].mvec.mvec[2], */
      /*        pointArray[j + i * NUM_POINT_IN_CIRC].mvec.mvec[3]); */
    }
  }
}

void rotateTorus(struct GaVector *pointArray, float angle3, float angle4) {
  struct GaBivector rotPlane3 = newGaBivec(5.0, 3.5, -4.0);
  struct GaRotor rotor3 = newGaRotor(angle3, rotPlane3);
  /* printf("rotor3:\n"); */
  /* printf("s: %f, e12: %f, e31: %f, e23: %f \n", rotor3.mvec.mvec[0], */
  /*        rotor3.mvec.mvec[4], rotor3.mvec.mvec[5], rotor3.mvec.mvec[6]); */
  struct GaBivector rotPlane4 = newGaBivec(-6.0, 2.0, 7.2);
  struct GaRotor rotor4 = newGaRotor(angle4, rotPlane4);
  /* printf("rotor4:\n"); */
  /* printf("s: %f, e12: %f, e31: %f, e23: %f \n", rotor4.mvec.mvec[0], */
  /*        rotor4.mvec.mvec[4], rotor4.mvec.mvec[5], rotor4.mvec.mvec[6]); */
  struct GaRotor rotor5 = {
      .mvec = GeometricProduct(rotor3.mvec, rotor4.mvec),
  };
  /* printf("rotor5:\n"); */
  /* printf("s: %f, e12: %f, e31: %f, e23: %f \n", rotor5.mvec.mvec[0], */
  /* rotor5.mvec.mvec[4], rotor5.mvec.mvec[5], rotor5.mvec.mvec[6]); */

  struct GaRotor rotor5Inv;
  rotor5Inv.mvec = Inverse(rotor5.mvec);

  for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {
    struct GaVector vecRot = {
        .mvec = GeometricProduct(
            GeometricProduct(rotor5Inv.mvec, pointArray[i].mvec), rotor5.mvec),
    };

    pointArray[i] = vecRot;
  }
}

void projectTorus(struct GaVector *pointArray) {

  // \[\vec{v}_\parallel=(\vec{v}\cdot\vec{B})\vec{B}^{-1} \]
  // plane of projecting
  struct GaBivector screenPlane = newGaBivec(0.0, 1.0, 0.0);
  struct GaBivector screenPlaneInv = {
      .mvec = Inverse(screenPlane.mvec),
  };

  for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {
    struct GaVector vectorProj = {
        .mvec =
            GeometricProduct(DotProduct(pointArray[i].mvec, screenPlane.mvec),
                             screenPlaneInv.mvec),
    };
  }
}

struct Pixel {
  uint8_t x;
  uint8_t y;
};

void displayTorus(struct GaVector *pointArray) {
  // Project into a rendering buffer
  // This buffer is an array of arrays
  // the array is a char array filled with " "
  // the outer array is the row of the donut
  // the inner array is the collum of the array
  // add "#" depentend on the projected donut
  struct Pixel renderBuffer[NUM_POINT_IN_TORUS];

  for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {

    renderBuffer[i].x = (int)ceilf(pointArray[i].mvec.mvec[1] * DONUT_WIDTH /
                                   (minorRadius * 2 + majorRadius)) +
                        DONUT_WIDTH + DONUT_WIDTH_BUFFER;
    renderBuffer[i].y =
        DONUT_HEIGTH - (int)ceilf(pointArray[i].mvec.mvec[3] * DONUT_HEIGTH /
                                  (minorRadius * 2 + majorRadius));
  }

  printf("\x1b[2J");
  for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {
    /* printf("x: %d, y: %d", renderBuffer[i].x, renderBuffer[i].y); */
    printf("\x1b[%d;%dH", renderBuffer[i].y, renderBuffer[i].x);
    printf("#");
  }
}

int main() {
  printf("\x1b[2J");
  printf("\x1b[?25l");
  printf("\x1b[0;0H");
  // Create torus
  struct GaVector pointArray[NUM_POINT_IN_TORUS]; // NUM_POINT_IN_TORIS
  makeTorus(pointArray);

  float angle3 = 0.0;
  float angle4 = 0.0;

  int num = 0;
  while (1) {

    struct GaVector torusArray[NUM_POINT_IN_CIRC]; // NUM_POINT_IN_TORIS
    memcpy(torusArray, pointArray, sizeof(torusArray));

    angle3 += TAU / 3;
    if (angle3 > TAU) {
      angle3 -= TAU;
    }
    /* printf("\x1b[0;0H"); */
    /* printf("%f\n", angle3); */

    angle4 += TAU / 4;
    if (angle4 > TAU) {
      angle4 -= TAU;
    }

    rotateTorus(pointArray, angle3, angle4);

    projectTorus(pointArray);

    displayTorus(pointArray);

    /* printf("\x1b[30;0H"); */
    /* num++; */
    /* printf("%d", num); */
  }
  return 0;
}
