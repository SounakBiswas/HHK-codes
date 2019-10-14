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
#define WARMUP 1000
#define MCSTEPS 1000
#define RELAX_STEPS 5
#define BINSIZE 10
#define SFACM 1
#define TAUMAX 20

int neigh[NSITES*6];
int sq_samples;
int sq_measure; // the current number of times sq have been measured.
int nsites,ncells;
int nautocorr;
int taumax;
double temp,beta;
double j2;
double sx[NSITES], sy[NSITES], sz[NSITES];
double sxqr[NSITES], syqr[NSITES], szqr[NSITES];
double sxqi[NSITES], syqi[NSITES], szqi[NSITES];
double energy;
double energy2;
double energy_bin, energy2_bin;
int lx,ly,lz;
int binsize;
char binfname[200],autofname[200],sfacnamepref[200],sfacname[200];
double sqasqbre[78*NCELLS];
double sqasqbim[78*NCELLS];
long long int nmeasure;
int binno;

