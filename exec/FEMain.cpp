#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "JacobiSolver.H"
#include <iostream>
#include <array>
using namespace std;

float sourceFunction(array<double, DIM> x)
{
  double val = -0.2;
  // region 1
  float Rsquared = (x[1] - 9) * (x[1] - 9) + x[0] * x[0];
  if (Rsquared > 25 && Rsquared < 36)
  {
    val = 1.5;
  }
  return val;
}

int main(int argc, char **argv)
{
  if (argc != 2)
  {
    cout << "this program takes one argument that is the .node and .ele file prefix\n";
    return 1;
  }
  string prefix(argv[1]);
  string nodeFile = prefix + ".node";
  string eleFile = prefix + ".ele";

  // --- load mesh ---
  FEGrid grid(nodeFile, eleFile);

  // --- build operator ---
  FEPoissonOperator<double> op(grid);

  // --- make a modifiable copy of the global stiffness matrix ---
  SparseMatrix<double> A = op.matrix();

  // (no symmetry check — matrix is now nonsymmetric by design)

  // --- compute source terms at centroids ---
  int nElements = grid.getNumElts();
  vector<double> sourceTerms(nElements);
  for (int i = 0; i < nElements; ++i)
  {
    auto centroid = grid.centroid(i);
    sourceTerms[i] = sourceFunction(centroid);
  }

  // --- build RHS (global) ---
  vector<double> rhs;
  op.makeRHS(rhs, sourceTerms);

  // --- prescribe non-zero Dirichlet BC via a lambda ---
  auto Phi_omega = [](const Node &n) -> double
  {
    // example boundary condition:
    return std::sin(M_PI * n.getPosition()[0]);
  };
  op.applyDirichletBC(A, rhs, Phi_omega);

  // --- solve the full system for phi at all nodes ---
  vector<double> phi;
  JacobiSolver<double> solver;
  double residual = solver.solve(phi, A, rhs, 1e-6, 1000);

  // --- write out VTK with boundary values baked in ---
  FEWrite(&grid, &phi, "solution.vtk");

  cout << "Final Solver residual was " << residual << endl;
  return 0;
}
