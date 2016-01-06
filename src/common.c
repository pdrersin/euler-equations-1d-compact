#include <stdio.h>
#include <math.h>
#include "common.h"

void ArrayToGnuplot( const char* fname, real* x, real *y, const int n) {
    int i;
    FILE* output;
    if( (output=fopen(fname,"w")) != NULL) {
	for(i=0;i<n;i++) {
	    fprintf(output,"%e %e\n",x[i],y[i]);
	}
	fclose(output);
    } else {
	fprintf(stderr,"ArrayToGnuplot: Error opening file `%s' for writing.\n",fname);
    }
}

void GetArrayMinMax( real* a, const int n, real* min, real* max) {
    int i;
    *min = -log(0);
    *max = log(0);
    for(i=0;i<n;i++) {
	*min = fmin(a[i],*min);
	*max = fmax(a[i],*max);
    }
}

