#include "common.h"
#include "compact.h"
#include "domain.h"
#include "deriv.h"
#include "rk.h"

/* Take a single timestep of the rk solver */
void euler(staggeredGrid* g) {
    int j,i;
    CalculateDerivs(g,g->U[0],g->R[0]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[0][i][j] = g->U[0][i][j]+g->dt*g->R[0][i][j];
	}
    }
    g->t += g->dt;
}

/* Take a single timestep of the rk solver */
void rk2(staggeredGrid* g) {
    int j,i;
    CalculateDerivs(g,g->U[0],g->R[0]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[1][i][j] = g->U[0][i][j]+g->dt*g->R[0][i][j];
	}
    }
    CalculateDerivs(g,g->U[1],g->R[1]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[0][i][j] = 0.5*g->U[0][i][j]+0.5*g->U[1][i][j]+0.5*g->dt*g->R[1][i][j];
	}
    }
    g->t += g->dt;
}

/* Take a single timestep of the rk solver */
void rk3(staggeredGrid* g) {
    int j,i;
    CalculateDerivs(g,g->U[0],g->R[0]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[1][i][j] = g->U[0][i][j]+g->dt*g->R[0][i][j];
	}
    }
    CalculateDerivs(g,g->U[1],g->R[1]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[2][i][j] = 0.75*g->U[0][i][j]+0.25*g->U[1][i][j]+0.25*g->dt*g->R[1][i][j];
	}
    }
    CalculateDerivs(g,g->U[2],g->R[2]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[0][i][j] = 1.0/3.0*g->U[0][i][j] + 2.0/3.0*(g->U[2][i][j] + g->dt*g->R[2][i][j]);
	}
    }
    g->t += g->dt;
}

/* Take a single timestep of the rk solver */
void rk4(staggeredGrid* g) {
    int j,i;
    CalculateDerivs(g,g->U[0],g->R[0]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[1][i][j] = g->U[0][i][j]+0.5*g->dt*g->R[0][i][j];
	}
    }
    CalculateDerivs(g,g->U[1],g->R[1]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[2][i][j] = g->U[0][i][j]+0.5*g->dt*g->R[1][i][j];
	}
    }
    CalculateDerivs(g,g->U[2],g->R[2]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[3][i][j] = g->U[0][i][j]+g->dt*g->R[2][i][j];
	}
    }
    CalculateDerivs(g,g->U[3],g->R[3]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[0][i][j] = 
		(-g->U[0][i][j] + g->U[1][i][j] +
		 2.0*g->U[2][i][j] + g->U[3][i][j] + 
		 0.5*g->dt*g->R[3][i][j])/3.0;
	}
    }
    g->t += g->dt;
}

/* Take a single timestep of the rk solver */
void rk45(staggeredGrid* g) {
    int j,i;
    CalculateDerivs(g,g->U[0],g->R[0]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[1][i][j] = g->U[0][i][j]+
		0.2*g->dt*g->R[0][i][j];
	}
    }
    CalculateDerivs(g,g->U[1],g->R[1]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[2][i][j] = g->U[0][i][j]+
		0.075*g->dt*g->R[0][i][j]+
		0.225*g->dt*g->R[1][i][j];
	}
    }
    CalculateDerivs(g,g->U[2],g->R[2]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[3][i][j] = g->U[0][i][j]+
		0.97777777777778*g->dt*g->R[0][i][j]-
		3.73333333333333*g->dt*g->R[1][i][j]+
		3.55555555555557*g->dt*g->R[2][i][j];
	}
    }
    CalculateDerivs(g,g->U[3],g->R[3]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[4][i][j] = g->U[0][i][j]+
		2.95259868922420*g->dt*g->R[0][i][j]-
		11.5957933241884*g->dt*g->R[1][i][j]+
		9.82289285169944*g->dt*g->R[2][i][j]-
		0.290809327846365*g->dt*g->R[3][i][j];
	}
    }
    CalculateDerivs(g,g->U[4],g->R[4]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[5][i][j] = g->U[0][i][j]+
		2.84627525252525*g->dt*g->R[0][i][j]-
		10.7575757575758*g->dt*g->R[1][i][j]+
		8.90642271774347*g->dt*g->R[2][i][j]+
		0.278409090909091*g->dt*g->R[3][i][j]-
		0.273531303602058*g->dt*g->R[4][i][j];
	}
    }
    CalculateDerivs(g,g->U[5],g->R[5]);
    for(i=0;i<g->nVar;i++) {
	for(j=0;j<g->nx;j++) {
	    g->U[0][i][j] = g->U[0][i][j]+
		0.0911458333333333*g->dt*g->R[0][i][j]+
		0.449236298292902*g->dt*g->R[2][i][j]+
		0.651041666666667*g->dt*g->R[3][i][j]-
		0.322376179245283*g->dt*g->R[4][i][j]+
		0.130952380952381*g->dt*g->R[5][i][j];
	}
    }
    g->t += g->dt;
}



