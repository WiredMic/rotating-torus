#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

// Include basic geometric algebra structs
#include "bivector.h"
#include "mvec.h"
#include "rotor.h"
#include "vector.h"

// define constants
#define TAU 6.283185
#define minorRadius 1.0
#define majorRadius 3.0
#define NUM_POINT_IN_CIRC 128
#define NUM_CIRC_IN_TORUS 128
#define NUM_POINT_IN_TORUS (NUM_POINT_IN_CIRC * NUM_CIRC_IN_TORUS)
#define DONUT_WIDTH 25
#define DONUT_WIDTH_BUFFER 10
#define DONUT_HEIGTH (DONUT_WIDTH / 2.0)
#define DONUT_HEIGTH_BUFFER 5
#define FPS 30

// pass by reference
// \[\vec{v}' = R_2^\dag(R_1^\dag \vec{v} R_1 + \vec{r})R_2\]
// \[ = (R_1R_2)^\dag \vec{v} R_1R_2 + R_2^\dag\vec{r}R_2\]
void makeTorus(struct GaVector *pointArray) {

  struct GaVector firstPoint = {
      .mvec.mvec[1] = minorRadius, // point
      .mvec.mvec[6] = 1,           // tangent plane
  };
  struct GaVector translate = newGaVec(majorRadius, 0, 0);

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
          .mvec.mvec[4] = vectorRot.mvec.mvec[4],
          .mvec.mvec[5] = vectorRot.mvec.mvec[5],
          .mvec.mvec[6] = vectorRot.mvec.mvec[6],
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
  struct GaMultivector tmpRotor5 = GeometricProduct(rotor3.mvec, rotor4.mvec);
  struct GaRotor rotor5 = {
      .mvec.mvec[0] = tmpRotor5.mvec[0] / Norm(tmpRotor5),
      .mvec.mvec[4] = tmpRotor5.mvec[4] / Norm(tmpRotor5),
      .mvec.mvec[5] = tmpRotor5.mvec[5] / Norm(tmpRotor5),
      .mvec.mvec[6] = tmpRotor5.mvec[6] / Norm(tmpRotor5),
  };
  // normalize the new rotor

  /* printf("rotor5:\n"); */
  /* printf("s: %f, e12: %f, e31: %f, e23: %f, length %f \n",
   * rotor5.mvec.mvec[0], */
  /*        rotor5.mvec.mvec[4], rotor5.mvec.mvec[5], rotor5.mvec.mvec[6], */
  /*        Norm(tmpRotor5)); */

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

struct Pixel {
  uint8_t x;
  uint8_t y;
  char brightness;
  float zBuffer;
};

void projectTorus(struct GaVector *pointArray, struct Pixel *renderBuffer) {

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

    /* printf("%d: e1: %f, e2: %f, e3: %f\n", i, vectorProj.mvec.mvec[1], */
    /*        vectorProj.mvec.mvec[2], vectorProj.mvec.mvec[3]); */

    renderBuffer[i].x = (int)ceilf(vectorProj.mvec.mvec[1] * DONUT_WIDTH /
                                   (minorRadius * 2 + majorRadius)) +
                        DONUT_WIDTH + DONUT_WIDTH_BUFFER;
    renderBuffer[i].y = DONUT_HEIGTH + DONUT_HEIGTH_BUFFER -
                        (int)ceilf(vectorProj.mvec.mvec[3] * DONUT_HEIGTH /
                                   (minorRadius * 2 + majorRadius));

    // In order to find the brigtness of the pixel
    // The way to find the brigtness is to make a unit light vector and then
    // take the wedge product of the unit tangent plane the a pixel. The volume
    // of the resulting trivector is the cosinus angle between the light vector
    // and the plane. 0 > then the surface is facing the light 0 < then the
    // surface is not facing Everything in between is some inbetween brigtness

    // tanget plane
    // the dual vector must go out

    // light vector
    struct GaVector lightRayDirection = newGaVec(1, 1, 1);
    float raylenght = Norm(lightRayDirection.mvec);
    struct GaVector lightRay = {
        .mvec.mvec[1] = lightRayDirection.mvec.mvec[1] / raylenght,
        .mvec.mvec[2] = lightRayDirection.mvec.mvec[2] / raylenght,
        .mvec.mvec[3] = lightRayDirection.mvec.mvec[3] / raylenght,
    };

    // There are 12 brightness values
    // The wedge product gives values form -1 to 1
    // ceilf((+ 1) * 6) gives decret values from 1 to 12;
    int brightness = (int)ceilf(
        (WedgeProduct(lightRay.mvec, pointArray[i].mvec).mvec[7] + 1) * 6);

    // .,-~:;=!*#$@
    switch (brightness) {
    case 1: // smallest
      renderBuffer[i].brightness = '.';
      break;
    case 2:
      renderBuffer[i].brightness = ',';
      break;
    case 3:
      renderBuffer[i].brightness = '-';
      break;
    case 4:
      renderBuffer[i].brightness = '~';
      break;
    case 5:
      renderBuffer[i].brightness = ':';
      break;
    case 6:
      renderBuffer[i].brightness = ';';
      break;
    case 7:
      renderBuffer[i].brightness = '=';
      break;
    case 8:
      renderBuffer[i].brightness = '!';
      break;
    case 9:
      renderBuffer[i].brightness = '*';
      break;
    case 10:
      renderBuffer[i].brightness = '#';
      break;
    case 11:
      renderBuffer[i].brightness = '$';
      break;
    case 12:
      renderBuffer[i].brightness = '@';
      break;
    default:
      /* renderBuffer[i].brightness = '#'; */
    }

    // z bufer
    // this is the rejected vector
    // \[\vec{v}_perp = (\vec{v}\wedge\vec{B})\vec{B}^{-1} \]

    struct GaVector vectorRejct = {
        .mvec =
            GeometricProduct(WedgeProduct(pointArray[i].mvec, screenPlane.mvec),
                             screenPlaneInv.mvec),
    };

    renderBuffer[i].zBuffer = vectorRejct.mvec.mvec[2];
  }
}

void displayTorus(struct Pixel *renderBuffer) {

  for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {
  }

  printf("\x1b[2J");
  for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {
    /* printf("i: %i,x: %d, y: %d\n", i, renderBuffer[i].x, renderBuffer[i].y);
     */
    printf("\x1b[%d;%dH", renderBuffer[i].y, renderBuffer[i].x);
    printf("%c\n", renderBuffer[i].brightness);
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

  while (1) {

    struct GaVector torusArray[NUM_POINT_IN_TORUS]; // NUM_POINT_IN_TORIS
    memcpy(torusArray, pointArray, sizeof(torusArray));
    for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {

      /* printf("%d: e1: %f, e2: %f, e3: %f\n", i, torusArray[i].mvec.mvec[1],
       */
      /*        torusArray[i].mvec.mvec[2], torusArray[i].mvec.mvec[3]); */
    }

    angle3 += TAU / (FPS * 512);
    if (angle3 > TAU) {
      angle3 -= TAU;
    }

    angle4 += TAU / (FPS * 258);
    if (angle4 > TAU) {
      angle4 -= TAU;
    }

    rotateTorus(pointArray, angle3, angle4);

    struct Pixel renderBuffer[NUM_POINT_IN_TORUS];

    projectTorus(torusArray, renderBuffer);

    displayTorus(renderBuffer);
    printf("\x1b[40;0H");
    printf("%f\n", angle3);
    printf("%f\n", angle4);
    usleep(1 / FPS * 1000000);
  }
  return 0;
}
