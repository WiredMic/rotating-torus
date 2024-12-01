// This line **must** come **before** including <time.h> in order to
// bring in the POSIX functions such as `clock_gettime() from <time.h>`!
#define _POSIX_C_SOURCE 199309L

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// Include basic geometric algebra structs
#include "bivector.h"
#include "mvec.h"
#include "rotor.h"
#include "vector.h"

// define function
/// Convert seconds to milliseconds
#define SEC_TO_MS(sec) ((sec) * 1000)
/// Convert nanoseconds to milliseconds
#define NS_TO_MS(ns) ((ns) / 1000000)

// define constants
#define TAU 6.283185
#define minorRadius 1.0
#define majorRadius 3.0
#define ANGULARV1 3.5
#define ANGULARV2 2.0
#define NUM_POINT_IN_CIRC 128
#define NUM_CIRC_IN_TORUS 128
#define NUM_POINT_IN_TORUS (NUM_POINT_IN_CIRC * NUM_CIRC_IN_TORUS)
// x direction
#define DONUT_WIDTH 25
#define DONUT_WIDTH_BUFFER 10
// y direction
#define DONUT_HEIGTH (DONUT_WIDTH / 2.0)
#define DONUT_HEIGTH_BUFFER 5
#define SCREEN_WIDTH 70
#define SCREEN_HEIGHT 30
#define FPS 2
#define DELTATIME

/// Get a time stamp in milliseconds.
uint64_t millis() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
  uint64_t ms = SEC_TO_MS((uint64_t)ts.tv_sec) + NS_TO_MS((uint64_t)ts.tv_nsec);
  return ms;
}

// pass by reference
// \[\vec{v}' = R_2^\dag(R_1^\dag \vec{v} R_1 + \vec{r})R_2\]
// \[ = (R_1R_2)^\dag \vec{v} R_1R_2 + R_2^\dag\vec{r}R_2\]
void makeTorus(struct GaVector *pointArray) {

  // First point all other point are rotated from
  struct GaVector firstPoint = {
      .mvec.mvec[1] = minorRadius, // point
      .mvec.mvec[6] = 1,           // tangent plane
  };

  // Translate vector to translate the circle from origin to the correct place
  // in torus
  struct GaVector translate = newGaVec(majorRadius, 0, 0);

  // angle and rotion plane for circle
  float angle1 = TAU / NUM_POINT_IN_CIRC;
  struct GaBivector rotPlane1 = newGaBivec(1.0, 0.0, 0.0);

  // angle and rotation plane for torus
  float angle2 = TAU / NUM_CIRC_IN_TORUS;
  struct GaBivector rotPlane2 = newGaBivec(0.0, 1.0, 0.0);

  // Iterate over all point in torus
  for (int i = 0; i < NUM_CIRC_IN_TORUS; i++) {

    // rotation around in the e3e1 plane
    struct GaRotor rotor2 = newGaRotor(angle2 * i, rotPlane2);

    // Rotate translation vector
    // R_2^\dag\vec{r}R_2
    struct GaVector translateRot = {
        .mvec = GeometricProduct(
            GeometricProduct(Reverse(rotor2.mvec), translate.mvec),
            rotor2.mvec),
    };

    for (int j = 0; j < NUM_POINT_IN_CIRC; j++) {
      // rotation around the e1e2 plane
      struct GaRotor rotor1 = newGaRotor(angle1 * j, rotPlane1);

      // The rotor that rotates from first point to the point in the torus
      // \[R_1R_2\]
      struct GaRotor rotor12 = {
          .mvec = GeometricProduct(rotor1.mvec, rotor2.mvec),
      };

      // rotate firste point + tangent plane
      // \[(R_1R_2)^\dag \vec{v} R_1R_2\]
      struct GaVector vectorRot = {
          .mvec = GeometricProduct(
              GeometricProduct(Reverse(rotor12.mvec), firstPoint.mvec),
              rotor12.mvec),
      };

      // translate the rotated point to the correct place in the torus
      // add point and tangent plane
      // \[ (R_1R_2)^\dag \vec{v} R_1R_2 + R_2^\dag\vec{r}R_2\]
      struct GaVector pointVec = {
          .mvec.mvec[1] = vectorRot.mvec.mvec[1] + translateRot.mvec.mvec[1],
          .mvec.mvec[2] = vectorRot.mvec.mvec[2] + translateRot.mvec.mvec[2],
          .mvec.mvec[3] = vectorRot.mvec.mvec[3] + translateRot.mvec.mvec[3],
          .mvec.mvec[4] = vectorRot.mvec.mvec[4],
          .mvec.mvec[5] = vectorRot.mvec.mvec[5],
          .mvec.mvec[6] = vectorRot.mvec.mvec[6],
      };

      // Add point + tangent plane to point array
      pointArray[j + i * NUM_POINT_IN_CIRC] = pointVec;
    }
  }
}

void rotateTorus(struct GaVector *pointArray, float angle3, float angle4) {

  struct GaBivector rotPlane3 = newGaBivec(5.0, 3.5, -4.0);
  struct GaRotor rotor3 = newGaRotor(angle3, rotPlane3);

  struct GaBivector rotPlane4 = newGaBivec(-6.0, 2.0, 7.2);
  struct GaRotor rotor4 = newGaRotor(angle4, rotPlane4);

  struct GaRotor rotor5 = {
      .mvec = GeometricProduct(rotor3.mvec, rotor4.mvec),
  };

  struct GaRotor rotor5Inv = {
      .mvec = Inverse(rotor5.mvec),
  };

  for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {
    struct GaVector vecRot = {
        .mvec = GeometricProduct(
            GeometricProduct(rotor5Inv.mvec, pointArray[i].mvec), rotor5.mvec),
    };

    pointArray[i] = vecRot;
  }
}

struct Screen {
  char brightness[SCREEN_HEIGHT][SCREEN_WIDTH]; // y size, x size
  float zbuffer[SCREEN_HEIGHT][SCREEN_WIDTH];   // y size, x size
};

struct Screen projectTorus(struct GaVector *pointArray) {

  // Make screen
  struct Screen screen;
  // make base string "     \0"
  char baseStr[SCREEN_WIDTH]; // x value
  for (int j = 0; j < SCREEN_WIDTH - 1; j++) {
    baseStr[j] = ' ';
  }
  baseStr[SCREEN_WIDTH - 1] = '\0'; // Scree width -1 ( because index from 0 )

  // copy base string to all strings in screen buffer
  for (int i = 0; i < SCREEN_HEIGHT; i++) { // y value
    strcpy(screen.brightness[i], baseStr);
  }

  for (int i = 0; i < SCREEN_HEIGHT; i++) {
    for (int j = 0; j < SCREEN_WIDTH; j++) {
      screen.zbuffer[i][j] = 0;
    }
  }

  // \[\vec{v}_\parallel=(\vec{v}\cdot\vec{B})\vec{B}^{-1} \]
  // plane of projecting
  struct GaBivector screenPlane = newGaBivec(0.0, 1.0, 0.0);
  struct GaBivector screenPlaneInv = {
      .mvec = Inverse(screenPlane.mvec),
  };

  // light vector
  struct GaVector lightRayDirection = newGaVec(1, 0.2, 1);
  float raylenght = Norm(lightRayDirection.mvec);
  struct GaVector lightRay = {
      .mvec.mvec[1] = lightRayDirection.mvec.mvec[1] / raylenght,
      .mvec.mvec[2] = lightRayDirection.mvec.mvec[2] / raylenght,
      .mvec.mvec[3] = lightRayDirection.mvec.mvec[3] / raylenght,
  };

  for (int i = 0; i < NUM_POINT_IN_TORUS; i++) {
    // Project point vector onto the screen plane
    struct GaVector vectorProj = {
        .mvec =
            GeometricProduct(DotProduct(pointArray[i].mvec, screenPlane.mvec),
                             screenPlaneInv.mvec),
    };

    // Scale the x and y value fit the sceeen
    int x = (int)ceilf(vectorProj.mvec.mvec[1] * DONUT_WIDTH /
                       (minorRadius * 2 + majorRadius)) +
            DONUT_WIDTH + DONUT_WIDTH_BUFFER;
    int y = DONUT_HEIGTH + DONUT_HEIGTH_BUFFER -
            (int)ceilf(vectorProj.mvec.mvec[3] * DONUT_HEIGTH /
                       (minorRadius * 2 + majorRadius));

    // In order to find the brigtness of the pixel
    // The way to find the brigtness is to make a unit light vector and then
    // take the wedge product of the unit tangent plane the a pixel. The volume
    // of the resulting trivector is the cosinus angle between the light vector
    // and the plane. 0 > then the surface is facing the light 0 < then the
    // surface is not facing Everything in between is some inbetween brigtness

    // The tanget plane was added at the start and has been rototed the same as
    // the point vectors

    // There are 12 brightness values
    // The wedge product gives values form -1 to 1
    // ceilf((+ 1) * 6) gives decret values from 1 to 12;
    int brightness = (int)ceilf(
        (WedgeProduct(lightRay.mvec, pointArray[i].mvec).mvec[7] + 1) * 6);

    char brightnessPixel;
    // .,-~:;=!*#$@
    switch (brightness) {
    case 1: // smallest
      brightnessPixel = '.';
      break;
    case 2:
      brightnessPixel = ',';
      break;
    case 3:
      brightnessPixel = '-';
      break;
    case 4:
      brightnessPixel = '~';
      break;
    case 5:
      brightnessPixel = ':';
      break;
    case 6:
      brightnessPixel = ';';
      break;
    case 7:
      brightnessPixel = '=';
      break;
    case 8:
      brightnessPixel = '!';
      break;
    case 9:
      brightnessPixel = '*';
      break;
    case 10:
      brightnessPixel = '#';
      break;
    case 11:
      brightnessPixel = '$';
      break;
    case 12:
      brightnessPixel = '@';
      break;
    default:
      brightnessPixel = ' ';
    }

    // z bufer
    // this is the rejected vector
    // \[\vec{v}_perp = (\vec{v}\wedge\vec{B})\vec{B}^{-1} \]
    struct GaVector vectorRejct = {
        .mvec =
            GeometricProduct(WedgeProduct(pointArray[i].mvec, screenPlane.mvec),
                             screenPlaneInv.mvec),
    };

    float zBuffer;
    if (vectorRejct.mvec.mvec[2] == 0) {
      zBuffer = vectorRejct.mvec.mvec[2] + 0.01;
    } else {
      zBuffer = vectorRejct.mvec.mvec[2];
    }

    if (zBuffer > screen.zbuffer[y][x] || screen.zbuffer[y][x] == 0) {

      screen.zbuffer[y][x] = zBuffer;

      screen.brightness[y][x] = brightnessPixel;
    }
  }

  return screen;
}

void displayTorus(struct Screen screen) {

  // place cursor at 1,1
  // clear screen
  // print strings

  printf("\x1b[1;1H"); // place cursor in 1 1
  printf("\x1b[2J");   // clear screen
  for (int i = 0; i < SCREEN_HEIGHT; i++) {
    printf("%s\n", screen.brightness[i]);
  }
}

int main() {
  // timing
  uint64_t startMillis = millis();

  printf("\x1b[?25l");
  // Create torus
  struct GaVector pointArray[NUM_POINT_IN_TORUS]; // NUM_POINT_IN_TORIS
  makeTorus(pointArray);

  float angle3 = 0.0;
  float angle4 = 0.0;

  while (1) {
    // time
    float deltat = ((float)millis() - (float)startMillis) / 1000; // seconds
    printf("\x1b[30;1H");
    printf("fps: %d\n", (int)(1 / deltat));
    startMillis = millis();

    struct GaVector torusArray[NUM_POINT_IN_TORUS]; // NUM_POINT_IN_TORIS
    memcpy(torusArray, pointArray, sizeof(torusArray));

    // add angle form angular velocity
    angle3 += ANGULARV1 * deltat;
    if (angle3 > TAU) {
      angle3 -= ((int)angle3 % (int)TAU) * TAU;
    }

    angle4 += ANGULARV2 * deltat;
    if (angle4 > TAU) {
      angle4 -= ((int)angle4 % (int)TAU) * TAU;
    }

    rotateTorus(torusArray, angle3, angle4);

    struct Screen screen = projectTorus(torusArray);

    displayTorus(screen);
  }
  return 0;
}
