#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "JacobiSolver.H"
#include <iostream>
#include <array>
using namespace std;

float sourceFunction(array<double, DIM> X)
{
  return X[0]*X[0]+X[1]*X[1];
}

int main(int argc, char **argv)
{
  // Parse command line arguments
  if(argc < 2 || argc > 3)
    {
      cout << "Usage:" << endl;
      cout << "  " << argv[0] << " <prefix> [max_area]" << endl;
      return 1;
    }

  string prefix(argv[1]);
  double max_area = 2.0;

  if (argc >= 2) {
    string arg2 = argv[2];
    max_area = stod(arg2);
  }

  // --- load mesh ---
  FEGrid grid;
  if (argc == 2)
  {
    string nodeFile = prefix+".node";
    string eleFile  = prefix+".ele";

    grid = FEGrid(nodeFile, eleFile);
  } else
  {
    string polyFile = prefix+".poly";
    
    grid = FEGrid(polyFile, max_area);
  }

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
    // return std::sin(M_PI * n.getPosition()[0]);
    // exact boundary condition
    return sourceFunction(n.getPosition());
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
