#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "ReInsert.H"
#include "JacobiSolver.H"
#include <iostream>
#include <array>
using namespace std;

float sourceFunction(array<double, DIM> x)
{
  double val=-.2;
  // region 1
  float Rsquared=(x[1]-9)*(x[1]-9)+x[0]*x[0];
  if(Rsquared > 25 && Rsquared < 36)
    {
      val = 1.5;
    }
  return val;
}

int main(int argc, char** argv)
{
  // Parse command line arguments
  if(argc < 2 || argc > 4)
    {
      cout << "Usage:" << endl;
      cout << "  " << argv[0] << " <prefix> [max_area] [--extrude]" << endl;
      return 1;
    }

  string prefix(argv[1]);
  bool extrude = false;
  double max_area = -1.0;

  if (argc >= 3) {
    string arg2 = argv[2];
    if (arg2 == "--extrude") {
      extrude = true;
    } else {
      max_area = stod(arg2);
    }
  }


  FEGrid grid;
  if (argc == 2 || (argc == 3 && extrude))
  {
    string nodeFile = prefix+".node";
    string eleFile  = prefix+".ele";

    grid = FEGrid(nodeFile, eleFile);
  } else
  {
    string polyFile = prefix+".poly";
    
    grid = FEGrid(polyFile, max_area);
  }

  FEPoissonOperator op(grid);

  const SparseMatrix& A = op.matrix(); 
  if(!A.symmetric())
    {
      cout<<" error in matrix assembly.  This should be a symmetric matrix\n";
      return 2;
    }
  array<int, 2> index = {{0,0}};
  cout << index[0] << " , " << index[1] << endl;
  
  for(unsigned int i=0; i<A.M(); i++, index[0]++, index[1]++)
    {
      if(A[index] <=0)
	{
	  cout <<"negative or zero diagonal detected "<<index[0]<<endl;
	}
    }

  int nElements = grid.getNumElts();
  vector<double> sourceTerms(nElements);
  array<double, DIM> centroid;
  for(int i=0; i<nElements; i++)
    {
      centroid = grid.centroid(i);
      sourceTerms[i] =sourceFunction(centroid);
    }

  vector<double> rhs, internalNodes, phi;
  op.makeRHS(rhs, sourceTerms);

  JacobiSolver solver;
  int iterations = 1000;
  double residual = solver.solve(internalNodes, op.matrix(), rhs, 1E-6, iterations);

  reinsert(grid, internalNodes, phi);

  FEWrite(&grid, &phi, "solution.vtk");

  cout<<" Final Solver residual was "<<residual<<endl;
  
  return 0;
  
}
