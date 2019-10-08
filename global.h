#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define LX 20
#define LY LX
#define LZ LX
#define NCELLS LX*LY*LZ 
#define NSITES 12*NCELLS
#define J2 1
#define TEMP 1.0
#define WARMUP 1000
#define MCSTEPS 1000
#define RELAX_STEPS 5

int neigh[NSITES*6];
int nsites,ncells;
double temp,beta;
double j2;
double sx[NSITES], sy[NSITES], sz[NSITES];
double energy;
double energy2;
int lx,ly,lz;

