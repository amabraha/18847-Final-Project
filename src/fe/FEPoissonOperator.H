#ifndef _FEPOISSONOPERATOR_H_
#define _FEPOISSONOPERATOR_H_
#include <cstdio> 
#include <cmath> 
#include <cassert>
#include <cmath>
#include <vector>
#include "Node.H"
#include "Element.H"
#include "FEGrid.H"
#include "SparseMatrix.H"
using namespace std;
/// Interface to the matrix representation of a finite element operator.
/**
This class constructs and stores a SparseMatrix object (m_matrix) using a finite element discretization defined by FEGrid.
 */
class FEPoissonOperator
{
 public:
  /// Construct from an FEGrid object.
  FEPoissonOperator(const FEGrid& a_grid);
  /// Construct a right-hand side, stored in a vector of length numInteriorNodes.
  void makeRHS(
               vector<double> & a_rhsAtNodes, 
               const vector<double> & a_FCentroids) const;

  /// const access to the FEGrid object.
  const FEGrid& getFEGrid() const;

  /// const access to the SparseMatrix object.
  const SparseMatrix& matrix() const;

 private:
  SparseMatrix m_matrix;
  FEGrid m_grid;
};
#endif
