#ifndef POISSONDISK_H
#define POISSONDISK_H
#include "matrix.h"
using namespace std;
typedef unsigned char uchar;

void PoissonDisk(int* Ir, int* Ig, int* Ib, vec3<int>* center, int numPixel, int numSample, float r);
void AutoPoissonDisk(int* Ir, int* Ig, int* Ib, vec3<int>* center, int numSpace, int numSample);

#endif
