#include <stdlib.h>
#include <stdio.h>
#include <rrlinalg.h>

#include "common.h"
#include "compact.h"


CompactWorkspace* AllocateCompactWorkspace( int n, int type, real* xm, real* x ) {
    int i;
    CompactWorkspace* cws = malloc(sizeof(CompactWorkspace));
    cws->n = n;
    cws->type = type;
    cws->xm = xm;
    cws->x = x;

    switch(type) {
	default:
	    fprintf(stderr,"Defaulting to Non-Periodic grid.\n");
	case NONPERIODIC:
	    fprintf(stderr,"Why don't you stick with a periodic grid for now.\n");
	    exit(1);
	break;
	case PERIODIC:
	    cws->result = malloc(sizeof(real)*n);
	    cws->work1 = malloc(sizeof(real)*n);
	    cws->work2 = malloc(sizeof(real)*n);
	    cws->A = malloc(sizeof(real*)*3);
	    for(i=0;i<3;i++) {
		cws->A[i] = malloc(sizeof(real)*n);
	    }
	break;
    }
    return cws;
}

void FreeCompactWorkspace( CompactWorkspace* cws) {
    int i;
    free(cws->result);
    if(cws->work1) free(cws->work1);
    if(cws->work2) free(cws->work2);
    for(i=0;i<3;i++) {
	free(cws->A[i]);
    }
    free(cws->A);
    free(cws);
}

void Compact6PeriodicDeriv( CompactWorkspace* cws, real* f ) {
    int i, ip32, ip12, im12, im32;
    int n = cws->n;
    real h = cws->xm[1] - cws->xm[0];

    for(i=0;i<n;i++) {
	cws->A[0][i] = 9.0/62.0;
	cws->A[1][i] = 1.0;
	cws->A[2][i] = 9.0/62.0;
	ip32 = (i+2+n)%(n);
	ip12 = (i+1  +n)%(n);
	im12 = (i-0+n)%(n);
	im32 = (i-1+n)%(n);
	cws->result[i] = 17.0/62.0*(f[ip32]-f[im32])/(3.0*h) + 63.0/62.0*(f[ip12]-f[im12])/h;
	//cws->result[i] = (f[ip12]-f[im12])/h;
    }

    rr_tridiag_periodic(cws->A, cws->result, n, cws->work1, cws->work2);
}

void Compact6PeriodicInterp( CompactWorkspace* cws, real* f ) {
    int i, ip32, ip12, im12, im32;
    int n = cws->n;

    for(i=0;i<n;i++) {
	cws->A[0][i] = 3.0/10.0;
	cws->A[1][i] = 1.0;
	cws->A[2][i] = 3.0/10.0;
	ip32 = (i+1+n)%(n);
	ip12 = (i  +n)%(n);
	im12 = (i-1+n)%(n);
	im32 = (i-2+n)%(n);
	cws->result[i] = 0.05*(f[ip32]+f[im32]) + 0.75*(f[ip12]+f[im12]);
    }
    rr_tridiag_periodic(cws->A, cws->result, n, cws->work1, cws->work2);
    //for(i=0;i<n;i++) {
	//fprintf(stderr,"f[%i] = %f,  interp[%i] = %f\n",i,f[i],i,cws->result[i]);
    //}
}

void Compact6PeriodicFilter( CompactWorkspace* cws, real* f ) {
    int i, ip2, ip1, im1, im2, ip3, im3, ip4, im4;
    int n = cws->n;

    real af, a0, a1, a2, a3, a4;
    af = 0.3;
    //a0 = 11.0/16.0 + 5.0*af/8.0;
    //a1 = 15.0/32.0 + 17.0*af/16.0;
    //a2 = -3.0/16.0 + 3.0*af/8.0;
    //a3 = 1.0/32.0 - af/16.0;
    a0 = (93.0+70.0*af)/128.0;
    a1 = (7.0+18.0*af)/16.0;
    a2 = (-7.0+14.0*af)/32.0;
    a3 = 1.0/16.0 - af/8.0;
    a4 = -1.0/128.0+af/64.0;
    //printf("a,b,c,d = %f %f %f %f\n",a,b,c,d);

    for(i=0;i<n;i++) {
	cws->A[0][i] = af;
	cws->A[1][i] = 1.0;
	cws->A[2][i] = af;
	ip4 = (i+4+n)%(n);
	ip3 = (i+3+n)%(n);
	ip2 = (i+2+n)%(n);
	ip1 = (i+1+n)%(n);
	im1 = (i-1+n)%(n);
	im2 = (i-2+n)%(n);
	im3 = (i-3+n)%(n);
	im4 = (i-4+n)%(n);
	cws->result[i] = a0*f[i] + 0.5*a1*(f[ip1]+f[im1]) + 0.5*a2*(f[ip2]+f[im2]) + 0.5*a3*(f[ip3]+f[im3]) + 0.5*a4*(f[ip4]+f[im4]);
    }
    rr_tridiag_periodic(cws->A, cws->result, n, cws->work1, cws->work2);
}
