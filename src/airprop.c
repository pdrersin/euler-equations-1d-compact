#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "airprop.h"

real WaveEqnC( real x ) {
    return 1.0;
}

// Assume constant R for now
real AirR() {
    return 287.0; // J/kg/K
}

// Sutherlands Formula for dynamic viscosity
real AirMu(real T) {
    return MU0_SUTH * (T0_SUTH+C_SUTH)/(T+C_SUTH) * (T/T0_SUTH) * sqrt(T/T0_SUTH);
}

real AirK(real T) {
    real T2 = T*T;
    return C0_K + C1_K*T + C2_K*T2 + C3_K*T2*T;
}

real AirCp(real T) {
    real ct = C_CP/T;
    real et = E_CP/T;
    real k1 = ct/sinh(ct);
    real k2 = et/sinh(et);
    return A_CP + B_CP*k1*k1 + D_CP*k2*k2;
}
real AirGam(real T) {
    //real Cp = AirCp(T);
    //real R = AirR();
    //return Cp / (Cp-R);
    return 1.4;
}


