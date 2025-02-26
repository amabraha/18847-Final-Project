
#ifndef _REINSERT_H_
#define _REINSERT_H_
#include "FEGrid.H"
#include <vector>

/// given a solution defined on the internal nodes of the grid, builds a full nodal solution with a value of 0 for the boundary nodes
void reinsert(const FEGrid& a_grid, const vector<double>& a_internalSolution, vector<double>& a_fullVector);

#endif
