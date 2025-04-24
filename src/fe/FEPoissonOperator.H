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
#include <array>
using namespace std;
/// Interface to the matrix representation of a finite element operator.
/**
This class constructs and stores a SparseMatrix object (m_matrix) using a finite element discretization defined by FEGrid.
 */
template <typename T>
class FEPoissonOperator
{
public:
  /// Construct from an FEGrid object.
  FEPoissonOperator(const FEGrid &a_grid);
  /// Construct a right-hand side, stored in a vector of length numInteriorNodes.
  void makeRHS(
      vector<T> &a_rhsAtNodes,
      const vector<T> &a_FCentroids) const;

  /// const access to the FEGrid object.
  const FEGrid &getFEGrid() const;

  /// const access to the SparseMatrix object.
  const SparseMatrix<T> &matrix() const;

  // apply non-zero Dirichlet boundary conditions
  void applyDirichletBC(SparseMatrix<T> &L,
                        std::vector<T> &rhs,
                        std::function<T(const Node &)> Phi_omega) const;

private:
  FEGrid m_grid;
  SparseMatrix<T> m_matrix;
};
#endif
