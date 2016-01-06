#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "common.h"
#include "initialconditions.h"
#include "airprop.h"

real rho0 = 1.225;
real p0 = 101325.0;
real k = 101325 * 0.01;
#define c sqrt(1.4*p0/rho0)

real P_Initial( real x ) {
    //return 101325.0 + 10000.0 * exp(-SQR((x-0.4)*20.0));
    return 101325.0 + k*WaveEqn_Initial(x);
}

real Rho_Initial( real x ) {
    return 1.225 + k/c/c * WaveEqn_Initial(x);
}

real U_Initial( real x ) {
    return 0.0 - k/c/rho0*WaveEqn_Initial(x);
}

real U1_Initial( real x ) {
    return Rho_Initial(x);
}

real U2_Initial( real x ) {
    return Rho_Initial(x) * U_Initial(x);
}

real U3_Initial( real x ) {
    real rho = Rho_Initial(x);
    real p   = P_Initial(x);
    real u   = U_Initial(x);
    real R   = AirR();
    return rho*u*u + p/(AirGam( p/rho/R )-1.0);
}

real WaveEqn_Initial( real x ) {
    //return 1.0 + exp(-SQR((x-0.2)*30.0));
    int i;
    real sum=0.0;
    srand48(4);
    for(i=2;i<201;i+=2) {
	sum += pow(i,-3.0/3.0)*sin(x*i*M_PI + drand48()*2.0*M_PI);
    }
    return sum;
}

void GetInitialConditions( real x, real* U ) {
#ifdef WAVE_EQN
    U[0] = WaveEqn_Initial( x );
#else
    U[0] = U1_Initial( x );
    U[1] = U2_Initial( x );
    U[2] = U3_Initial( x );
#endif
}

