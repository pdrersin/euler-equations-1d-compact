#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "common.h"
#include "compact.h"
#include "domain.h"
#include "airprop.h"
#include "deriv.h"

void CalculateCellCenteredFluxes( real** F, real** U, const int nx ) {
    int i;
    real gam, T, p;
#ifdef WAVE_EQN
    for(i=0;i<nx;i++) {
	F[0][i] = -WaveEqnC(0)*U[0][i];
    }
#else
    for(i=0;i<nx;i++) {
	// APPROXIMATE: gamma=1.4
	p = (U[2][i] - SQR(U[1][i])/U[0][i])*0.4;
	T = p/U[0][i]/AirR();
	gam = AirGam(T);
	F[0][i] = U[1][i];
	F[1][i] = (2.0-gam)*SQR(U[1][i])/U[0][i] + (gam-1.0)*U[2][i];
	F[2][i] = gam*U[1][i]*U[2][i]/U[0][i] + (1.0-gam)*CUBE(U[1][i])/SQR(U[0][i]);
    }
#endif
}

void CalculateDerivs(staggeredGrid* g, real** input, real** deriv) {
    int i;
    // Cell-centered fluxes (g->F) FROM cell-centered state vector (input)
    CalculateCellCenteredFluxes( g->F, input, g->nx );

    //Loop through all variables and calculate dF/dx at cell-centers
    for(i=0;i<g->nVar;i++) {
	// Interpolate those fluxes on to cell-faces
	Compact6PeriodicInterp( g->compact, g->F[i] );
	memcpy(g->F[i], g->compact->result, sizeof(real)*g->nx);

	// Differentiate with staggered scheme to return to cell-centers
	Compact6PeriodicDeriv( g->compact, g->F[i] );
	memcpy(deriv[i], g->compact->result, sizeof(real)*g->nx);
    }
}

void FilterAll(staggeredGrid* g) {
    int i;
    for(i=0; i<g->nVar;i++) {
	Compact6PeriodicFilter( g->compact, g->U[0][i] );
	memcpy(g->U[0][i], g->compact->result, sizeof(real)*g->nx);
    }
}
