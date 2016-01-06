#ifndef __INITIALCONDITIONS_H__
#define __INITIALCONDITIONS_H__


real WaveEqn_Initial( real x );

// Initial primitive variables
real P_Initial( real x );
real Rho_Initial( real x );
real U_Initial( real x );


// Initial conserved variables
real U1_Initial( real x );
real U2_Initial( real x );
real U3_Initial( real x );

// Return a full vector of initial conditions
void GetInitialConditions( real x, real* U );

#endif /*__INITIALCONDITIONS_H__*/
