#include <stdlib.h>
#include <stdio.h>
#define LX 4
#define LY LX
#define LZ LX
#define NCELLS LX*LY*LZ 
#define NSITES 12*NCELLS
#define J2 1
#define TEMP 1.0

int neigh[NSITES*6];
int nsites,ncells;
double temp,beta;
double j2;
double sx[NSITES], sy[NSITES], sz[NSITES];
int lx,ly,lz;

