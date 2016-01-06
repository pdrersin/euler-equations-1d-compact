#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "common.h"
#include "glinterface.h"
#include "gldisplay.h"
#include "airprop.h"
#include "compact.h"
#include "domain.h"
#include "rk.h"
#include "main.h"

staggeredGrid* grid;


int main(int argc, char **argv) {

    int nx = 201;
    int nVar = 3;
    int nLevels = 6;
    int nDerivs = 6;
    
    grid = AllocStaggeredGrid( nx, nVar, nLevels, nDerivs );
    grid->RK = rk45;

    
    RegularEdgeCoordinates( grid, 0.0, 1.0 );
    //ReadEdgeCoordinates( grid, "../dat/points.dat" );
    InitializeStaggeredGrid( grid );
    grid->dt = 0.00001;
    CalculatePrimitiveVariables( grid );

    update_plan();
    update_fft();

    //OutputStaggeredGrid( grid, "../dat/grid.dat", GNUPLOT );

    initGlut( argc, argv );
    glutMainLoop();


    return 0;
}

