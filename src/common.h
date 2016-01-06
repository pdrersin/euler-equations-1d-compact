#ifndef __COMMON_H__
#define __COMMON_H__

typedef double real;

//#define WAVE_EQN

#ifndef MIN
#define MIN(A,B) ((A)<(B)?(A):(B))
#endif

#ifndef MAX
#define MAX(A,B) ((A)<(B)?(B):(A))
#endif

#define BOUND(A,B,C) MIN(MAX(A,B),C)
#define ABS(A) ((A)<0?(-(A)):(A))
#define SQR(A) ((A)*(A))
#define CUBE(A) ((A)*(A)*(A))
#define SIGN(A) ((A<0.0)?(-1):(1))

void ArrayToGnuplot( const char* fname, real* x, real* y, const int n );

void GetArrayMinMax( real* a, const int n, real* min, real* max);

#endif /*__COMMON_H__*/
