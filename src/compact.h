#ifndef __COMPACT_H__
#define __COMPACT_H__

#define NONPERIODIC 0
#define PERIODIC 1

typedef struct {
    int type;            // Periodic vs. Non-Periodic
    int n;               // # grid points
    real *xm, *x;        // Cell-centered and face-centered grids
    real** A;
    real *result;
    real *work1, *work2;
} CompactWorkspace;

CompactWorkspace* AllocateCompactWorkspace( int n, int type, real* xm, real* x );
void FreeCompactWorkspace( CompactWorkspace* cws);
void Compact6PeriodicDeriv( CompactWorkspace* cws, real* f );
void Compact6PeriodicInterp( CompactWorkspace* cws, real* f );
void Compact6PeriodicFilter( CompactWorkspace* cws, real* f );

#endif /*__COMPACT_H__*/
