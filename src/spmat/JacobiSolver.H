#ifndef _JACOBISOLVER_H_
#define _JACOBISOLVER_H_

#include "SparseMatrix.H"
#include <vector>

class JacobiSolver
{
public:
  /// Solves the linear system a_A*a_phi = a_rhs  using Point Jacobi iteration.
  /**
  Point Jacobi iteration is applied for a_iter iterations, or until 
  \f$ max_j \f$  |(a_A* a_phi)[j] - a_rhs[j]| \f$ \leq \f$  a_tolerance *  \f$ max_j \f$ |a_rhs[j]| , whichever comes first.
   */
  double solve(
               vector<double>& a_phi,
               const SparseMatrix& a_A, 
               const vector<double>& a_rhs, 
               const double& a_tolerance, 
               int a_iter);
};

#endif
