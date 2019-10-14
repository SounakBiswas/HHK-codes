#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#define LX 4
#define LY LX
#define LZ LX
#define NCELLS LX*LY*LZ 
#define NSITES 12*NCELLS
#define J2 1
#define TEMP 1.0
#define WARMUP 10000
#define MCSTEPS 10000
#define RELAX_STEPS 5
#define BINSIZE 10
#define TAUMAX 20

int neigh[NSITES*6];
int nsites,ncells;
int nautocorr;
int taumax;
double temp,beta;
double j2;
double sx[NSITES], sy[NSITES], sz[NSITES];
double energy;
double energy2;
double energy_bin, energy2_bin;
int lx,ly,lz;
int binsize;
char binfname[200],autofname[200],fourierfname[200];
long long int nmeasure;
int binno;

