#ifndef __DERIV_H__
#define __DERIV_H__

void CalculateCellCenteredFluxes( real** F, real** U, const int nx );

void CalculateDerivs(staggeredGrid* g, real** input, real** deriv);

void FilterAll(staggeredGrid* g);

#endif /*__DERIV_H__*/
