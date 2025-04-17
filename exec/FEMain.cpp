#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "ReInsert.H"
#include "JacobiSolver.H"
#include <iostream>
#include <array>
using namespace std;

float sourceFunction(array<double, DIM> x)
{
  double val = -.2;
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
    cout << "this program takes one argument that is the .node and .ele ";
    cout << "file prefix" << endl;
    return 1;
  }
  string prefix(argv[1]);
  string nodeFile = prefix + ".node";
  string eleFile = prefix + ".ele";

  FEGrid grid(nodeFile, eleFile);

  FEPoissonOperator<double> op(grid);

  const SparseMatrix<double> &A = op.matrix();
  if (!A.symmetric())
  {
    cout << " error in matrix assembly.  This should be a symmetric matrix\n";
    return 2;
  }
  array<int, 2> index = {{0, 0}};
  cout << index[0] << " , " << index[1] << endl;

  for (unsigned int i = 0; i < A.M(); i++, index[0]++, index[1]++)
  {
    if (A[index] <= 0)
    {
      cout << "negative or zero diagonal detected " << index[0] << endl;
    }
  }

  int nElements = grid.getNumElts();
  vector<double> sourceTerms(nElements);
  array<double, DIM> centroid;
  for (int i = 0; i < nElements; i++)
  {
    centroid = grid.centroid(i);
    sourceTerms[i] = sourceFunction(centroid);
  }

  vector<double> rhs, internalNodes, phi;
  op.makeRHS(rhs, sourceTerms);

  JacobiSolver<double> solver;
  int iterations = 1000;
  double residual = solver.solve(internalNodes, op.matrix(), rhs, 1E-6, iterations);

  reinsert(grid, internalNodes, phi);

  FEWrite(&grid, &phi, "solution.vtk");

  cout << " Final Solver residual was " << residual << endl;

  return 0;
}
