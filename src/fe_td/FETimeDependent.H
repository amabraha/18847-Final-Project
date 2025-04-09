#ifndef _FETIMEDEPENDENT_H_
#define _FETIMEDEPENDENT_H_
#include <cstdio> 
#include <cmath> 
#include <cassert>
#include <cmath>
#include <vector>
#include <functional>
#include "Node.H"
#include "Element.H"
#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "SparseMatrix.H"
#include "JacobiSolver.H"
#include "ReInsert.H"
using namespace std;

/**
This class will perform the time dependent FEM solving for the equation dPhi/dt = L(Phi)-f
It should be provided both L and f and also intial and boundary conditions for Phi
 */
class FETimeDependent
{
 public:
  /// Construct solver from sparse matrix L and f
  /// f needs to be a function that takes in a time and returns a vector because it is time dependent
  FETimeDependent(const SparseMatrix& a_L, const function<vector<double>(double)>& a_f, const FEGrid& a_grid);

  /// solves for phi(time) using time steps of <= dt
  /// will populate a_phi_out with the output.
  void solve(double time, double dt, vector<double>& a_phi_out, 
             const vector<double>& a_initial_cond, const function<vector<double>(double)>& a_boundary_cond);

  /// solves for phi(time) using time steps of <= dt
  /// will populate a_phi_out with the output.
  /// will also writes sequence vtk data to a_filenamei.vtk
  void solve_write(double time, double dt, vector<double>& a_phi_out, 
                   const vector<double>& a_initial_cond, const function<vector<double>(double)>& a_boundary_cond,
                   string a_filename);

  


 private:
  SparseMatrix m_L;
  function<vector<double>(double)> m_f;
  FEGrid m_grid;

  ///steps a_phi_out from time to time+dt
  void step(double time, double dt, vector<double>& a_phi_out);
};
#endif
