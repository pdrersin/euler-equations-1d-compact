#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "array.h"
#include "initialconditions.h"
#include "airprop.h"
#include "compact.h"
#include "domain.h"

// Keep it simple and just allocate the plain grid
staggeredGrid* AllocStaggeredGrid( const int nx, const int nVar, const int nLevels, const int nDerivs) {
    staggeredGrid* newGrid;
    newGrid = malloc(sizeof(staggeredGrid));
    
    // Save array sizes
    newGrid->nx = nx;
    newGrid->nVar = nVar;
    newGrid->nLevels = nLevels;
    newGrid->nDerivs = nDerivs;

    // Allocate coordinate arrays
    newGrid->x  = (real*)Alloc1DArray( sizeof(real), nx+1);
    newGrid->xm = (real*)Alloc1DArray( sizeof(real), nx  );

    // Allocate conserved variables   
    newGrid->W = (real**)Alloc2DArray( sizeof(real), nVar, nx );
    newGrid->U = (real***)Alloc3DArray( sizeof(real), nLevels, nVar, nx );
    newGrid->R = (real***)Alloc3DArray( sizeof(real), nDerivs, nVar, nx );
    newGrid->F = (real**)Alloc2DArray( sizeof(real), nVar, nx );

    // Allocate the compact workspace
    newGrid->compact = AllocateCompactWorkspace( nx, PERIODIC, newGrid->xm, newGrid->x );

    // Return pointer to the grid
    return newGrid;
}

void FreeStaggeredGrid( staggeredGrid* grid ) {

    // Free coordinate arrays
    Free1DArray( (void*)grid->x  );
    Free1DArray( (void*)grid->xm );

    // Free conserved variable arrays
    Free2DArray( (void**)grid->W, grid->nVar );
    Free3DArray( (void***)grid->U, grid->nLevels, grid->nVar );
    Free3DArray( (void***)grid->R, grid->nDerivs, grid->nVar );
    Free2DArray( (void*)grid->F, grid->nVar );

    // Free compact scheme workspace
    FreeCompactWorkspace( grid->compact );

    // Free the whole grid
    free( grid );
}

void RegularEdgeCoordinates( staggeredGrid* grid, real xmin, real xmax) {
    int i;
    for(i=0;i<grid->nx+1;i++) {
	grid->x[i] = xmin + (xmax-xmin)*i/grid->nx;
    }

    for(i=0;i<grid->nx+1;i++) {
	grid->xm[i] = 0.5*(grid->x[i]+grid->x[i+1]);
    }
}
void ReadEdgeCoordinates( staggeredGrid* grid, const char* inputFileName ) {
    int i;
    int nPoints=0;
    char buffer[128];
    double din1;

    FILE* input;

    if( (input=fopen(inputFileName,"r")) != NULL ) {

	while( fgets(buffer,sizeof(buffer),input) ) {
	    if(buffer[0] != '#') {
		sscanf( buffer, "%le", &din1 );
		if( nPoints < grid->nx+1 ) {
		    grid->x[nPoints] = (real)din1;
		} else {
		    fprintf(stderr,"ReadEdgeCoordinates: Warning: Ignoring extra coordinate data.\n");
		    break;
		}
		nPoints++;
	    }
	}

	for(i=0;i<grid->nx+1;i++) {
	    grid->xm[i] = 0.5*(grid->x[i]+grid->x[i+1]);
	}

    } else {
	fprintf( stderr,"ReadEdgeCoordinates: Error opening file `%s' for reading.\n", inputFileName );
	exit(1);
    }
}

void InitializeStaggeredGrid( staggeredGrid* grid ) {
    int i,j;
    real* U0 = malloc(sizeof(real)*grid->nVar);

    // Loop through cell centers
    for(i=0;i<grid->nx;i++) {

	// Get a vector of initial conditions in the cell
	GetInitialConditions( grid->xm[i], U0 );

	for(j=0;j<grid->nVar;j++) {
	    grid->U[0][j][i] = U0[j];
	}
    }

    free( U0 );
}

void OutputStaggeredGrid( staggeredGrid* grid, const char* outputFileName, const int outputType ) {
    int i;
    FILE* output;
#ifndef WAVE_EQN
    int j;
    real U1, U2, U3, rho, u, p, T;
#endif

    if( (output=fopen(outputFileName,"w")) != NULL) {

#ifdef WAVE_EQN
	if( outputType == TECPLOT) {
	    fprintf( output, "VARIABLES = \"x\", \"y\"");
	    fprintf( output, "ZONE T=Grid, I=%i\n",grid->nx );
	} else {
	    fprintf( output, "# x y\n");
	}

	for(i=0;i<grid->nx;i++) {
	    fprintf(output,"%le %le\n", grid->xm[i], grid->U[0][0][i]);
	}
#else
	if( outputType == TECPLOT) {
	    fprintf( output, "VARIABLES = \"x\", \"rho\", \"u\", \"p\"");
	    fprintf( output, "ZONE T=Grid, I=%i\n",grid->nx );
	} else {
	    fprintf( output, "# x rho u p\n");
	}

	for(i=0;i<grid->nx;i++) {

	    // Conserved variables
	    U1 = grid->U[0][0][i];
	    U2 = grid->U[0][1][i];
	    U3 = grid->U[0][2][i];

	    // Primitive variables
	    rho = U1;
	    u = U2/U1;
	    
	    //   UH-OH ...  Let's just guess gamma = 1.4
	    p = ( 0.4 ) * (U3 - U2*U2/U1);
	    
	    // Fixed-point iteration to improve the estimate of p based on the correct gamma...
	    for(j=0;j<3;j++) {
		T = p/rho/AirR();
		p = (AirGam(T)-1.0)*(U3-U2*U2/U1);
	    }

	    fprintf(output,"%le %le %le %le\n", grid->xm[i], rho, u, p);
	}
#endif

    } else {
	fprintf(stderr,"OutputStaggerdGridToGnuplot: Error: Could not open `%s' for writing.\n",outputFileName);
    }
}

real GetCFLNumber( staggeredGrid* grid ) {
    int i;
    real u,c,p,T,gam, max=0.0;
    real** U = grid->U[0];
    for(i=0;i<grid->nx;i++) {
	u = U[1][i]/U[0][i];
	p = (U[2][i] - SQR(U[1][i])/U[0][i])*0.4;
	T = p/U[0][i]/AirR();
	gam = AirGam(T);
	//p = (U[2][i] - SQR(U[1][i])/U[0][i])*(gam-1.0);
	//T = p/U[0][i]/AirR();
	//gam = AirGam(T);
	c = sqrt(gam*AirR()*T);
	max = fmax(fabs(u+c),max);
	max = fmax(fabs(u-c),max);
    }
    return max * grid->dt / (grid->x[1]-grid->x[0]);
}

void CalculatePrimitiveVariables( staggeredGrid* grid ) {
    int i;
    real** U = grid->U[0];
    for(i=0; i<grid->nx; i++) {
	// rho ( just copy ):
	grid->W[0][i] = U[0][i];
	// u:
	grid->W[1][i] = U[1][i]/U[0][i];
	// p:
	grid->W[2][i] = (U[2][i] - SQR(U[1][i])/U[0][i])*0.4;
	//grid->W[0][i] = U[0][i];
	//grid->W[1][i] = U[1][i];
	//grid->W[2][i] = U[2][i];
    }
}


