#ifndef __AIRPROP_H__
#define __AIRPROP_H__

#define MW_AIR (28.97)

// Viscosity parameters 
// Wikipedia...  ouch.
#define C_SUTH (120.0)
#define T0_SUTH (291.15)
#define MU0_SUTH (18.27e-6)

// Thermal conductivity parameters
// Not sure about this, but it's from
// http://users.wpi.edu/~ierardi/PDF/air_k_plot.pdf ...
#define C0_K (-3.9333e-4)
#define C1_K (1.0184e-4)
#define C2_K (-4.8574e-8)
#define C3_K (1.5207e-11)

// Specific Heat at Constant Pressure parameters (Kuo, p. 641)
#define A_CP (9.9959e2)
#define B_CP (3.2413e2)
#define C_CP (3.0120e3)
#define D_CP (2.6165e2)
#define E_CP (1.4840e3)


real WaveEqnC( real x );

real AirR();
real AirMu(real T);
real AirK(real T);
real AirGam(real T);
real AirCp(real T);


#endif /*__AIRPROP_H__*/
