#include "JacobiSolver.H"
#include <cassert>
#include <iostream>
#include <array>
#include <complex>
using namespace std;

template <typename T>
double norm(const vector<T> &a_v)
{
  double rtn = 0;
  for (unsigned int i = 0; i < a_v.size(); i++)
  {
    double x = abs(a_v[i]);
    if (rtn < x)
      rtn = x;
  }
  return rtn;
}

// recursive helper for Jacobi point iteration
template <typename T>
double Jacobi_helper(
    vector<T> &a_phi,
    const SparseMatrix<T> &a_A,
    const vector<T> &a_rhs,
    const double &a_tolerance,
    double lambda,
    int a_curr_iter,
    int a_iter)
{
  vector<T> residual = vector<T>(a_phi.size(), T(0.0));
  // result of multiplying a_A with a_phi
  vector<T> La = a_A * a_phi;

  for (int indx = 0; indx < a_rhs.size(); indx++)
  {
    residual[indx] = a_rhs[indx] - La[indx];
  }

  // if we have converged or reached max number of iteratoins
  if (a_curr_iter == a_iter || norm(residual) <= a_tolerance * norm(a_rhs))
  {
    /*printf("initial norm(rhs): %f\nfinal norm(residual)/norm(rhs): %e\nnumber iterations: %d\n",
           norm(a_rhs),
           norm(residual) / norm(a_rhs),
           a_curr_iter);*/
    return norm(residual) / norm(a_rhs);
  }

  // else update a_phi for next iteration
  for (int indx = 0; indx < a_phi.size(); indx++)
  {
    a_phi[indx] = a_phi[indx] + T(lambda) * residual[indx];
  }

  return Jacobi_helper(
      a_phi,
      a_A,
      a_rhs,
      a_tolerance,
      lambda,
      a_curr_iter + 1,
      a_iter);
}

template <typename T>
double JacobiSolver<T>::solve(
    vector<T> &a_phi,
    const SparseMatrix<T> &a_A,
    const vector<T> &a_rhs,
    const double &a_tolerance,
    int a_iter)
{
  // check that A is square and a_rhs has same dimension
  assert(a_A.M() == a_A.N());
  assert(a_A.N() == a_rhs.size());

  // resize a_phi to have same dimension
  a_phi.resize(a_A.N());

  // zero entries in a_phi
  for (int indx = 0; indx < a_phi.size(); indx++)
  {
    a_phi[indx] = T(0.0);
  }

  // arbitrarily chosen from writeup
  double relaxation_parameter = 0.85;

  // compute max diagonal element in input matrix
  // note that at least one diagonal element must be positive for the matrix to have all positive eigenvalues
  double max_Lkk = 0;
  for (int k = 0; k < a_A.N() && k < a_A.M(); k++)
  {
    array<int, 2> diag_indx = array<int, 2>{k, k};
    double abs_val = abs(a_A[diag_indx]);
    if (max_Lkk < abs_val)
      max_Lkk = abs_val;
  }

  // call recursive helper
  return Jacobi_helper(
      a_phi,
      a_A,
      a_rhs,
      a_tolerance,
      relaxation_parameter / max_Lkk,
      1,
      a_iter);
}

// Explicit template instantiations
template class JacobiSolver<float>;
template class JacobiSolver<double>;
template class JacobiSolver<std::complex<float>>;
template class JacobiSolver<std::complex<double>>;
