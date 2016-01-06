#ifndef __DOMAIN_H__
#define __DOMAIN_H__

// Data structure to store a 1-dimensional staggered grid
// with staggered data arranged as per Nagarajan, Lele,
// and Ferziger, "A robust high-order compact method for 
// large eddy simulation."
//
// Conserved quantities are stored at cell centers, fluxes
// at the cell boundaries.  So conserved quantities are:
//
// |     |     |     |     |     |     |     |
// |  0  |  1  |  2  | ... | nx-2| nx-1| nx  |
// |     |     |     |     |     |     |     |
// 
// While the fluxes are:
//
// |     |     |     |     |     |     |     |
// 0     1     2     3    ...   nx-1  nx    nx+1 
// |     |     |     |     |     |     |     |
//
// No ghost cells are used.  Insted, only the compact schemes
// are altered in the neighborhood of the boundaries.  Of
// course the whole point is that the grid is not uniform, so
// schemes must be constructed on the fly.  `x' stores the
// edge coordinates, while xm stores the cell-centerd coords.

#define TECPLOT 0
#define GNUPLOT 1

typedef struct staggeredGrid_s {
    real **W;
    real ***U;
    real ***R;
    real **F;
    real *x;
    real *xm;
    real dt;
    real t;
    int nx;              // Number of CELLS in domain
    int nVar;            // Length of state vector
    int nLevels;         // Fractional step levels
    int nDerivs;         // Number of derivatives we need to store for RK
    CompactWorkspace* compact; // Avoid reallocating every time we need a derivative
    void (*RK)(struct staggeredGrid_s*);   // Pointer to RK solver function
} staggeredGrid;

// Keep it simple and just allocate the plain grid
staggeredGrid* AllocStaggeredGrid( const int nx, const int nVar, const int nLevels, const int nDerivs);

void FreeStaggeredGrid( staggeredGrid* grid );

void ReadEdgeCoordinates( staggeredGrid* grid, const char* inputFileName );

void RegularEdgeCoordinates( staggeredGrid* grid, real xmin, real xmax);

void InitializeStaggeredGrid( staggeredGrid* grid );

void OutputStaggeredGrid( staggeredGrid* grid, const char* outputFileName, const int outputType );

real GetCFLNumber( staggeredGrid* grid );

void CalculatePrimitiveVariables( staggeredGrid* grid );

#endif /*__DOMAIN_H__*/
