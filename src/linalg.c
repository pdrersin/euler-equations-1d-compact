#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "tdiag.h"

void tdiag(real** A, real* b, const int n) {
    int i=0;
    real fac;
    //Eliminate
    for(i=n-2; i>=0; i--) {
	fac = A[i][0]/A[i+1][1];
	A[i][1] = A[i][1]-fac*A[i+1][0];
	b[i] = b[i]-fac*b[i+1];
    }
    //Back-substitute
    b[0]=b[0]/A[0][1];
    for(i=1;i<n;i++) {
	b[i] = (b[i]-A[i][0]*b[i-1])/A[i][1];
    }
}

void pentadiag(real** A, real* b, const int n) {
    int i;
    real fac;
    //Eliminate 1
    for(i=2;i<n;i++) {
	fac=A[i][0]/A[i-1][1];
	if(A[i-1][1]==0) { fprintf(stderr,"Pentadiag failed due to zero pivot element.  May I suggest you get out and walk?\n"); exit(1); }
	A[i][1] -= fac*A[i-1][2];
	A[i][2] -= fac*A[i-1][3];
	A[i][3] -= fac*A[i-1][4];
	b[i] -= fac*b[i-1];
    }
    //Eliminate 2
    for(i=1;i<n;i++) {
	fac=A[i][1]/A[i-1][2];
	if(A[i-1][2]==0) { fprintf(stderr,"Pentadiag failed on elimination #2.  It's gonna be a long night.\n"); exit(1); }
	A[i][2] -= fac*A[i-1][3];
	A[i][3] -= fac*A[i-1][4];
	b[i] -= fac*b[i-1];
    }

    //Back-substitute
    b[n-1]/=A[n-1][2];
    b[n-2]=(b[n-2]-A[n-1][3]*b[n-1])/A[n-1][2];
    for(i=n-3;i>=0;i--) {
	b[i] = (b[i]-A[i][3]*b[i+1]-A[i][4]*b[i+2])/A[i][2];
    }
}
