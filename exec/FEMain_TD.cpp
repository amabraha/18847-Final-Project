#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "FETimeDependent.H"
#include "ReInsert.H"
#include "JacobiSolver.H"
#include <iostream>
#include <array>
using namespace std;

float sourceFunction(double time, array<double, DIM> x)
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
  if(argc != 2)
    {
      cout << "this program takes one argument that is the .node and .ele ";
      cout << "file prefix"<<endl;
      return 1;
    }
  string prefix(argv[1]);
  string nodeFile = prefix+".node";
  string eleFile  = prefix+".ele";

  FEGrid grid(nodeFile, eleFile);

  FEPoissonOperator op(grid);
  

  vector<double> initial_conditions(grid.getNumInteriorNodes());

  function<vector<double>(double)> rhs_f = [&grid, &op](double time)
  {

    int nElements = grid.getNumElts();
    vector<double> sourceTerms(nElements);

    array<double, DIM> centroid;

    for(int i=0; i<nElements; i++)
    {
    centroid = grid.centroid(i);
    sourceTerms[i] = sourceFunction(time, centroid);
    }
  
    vector<double> rhs;
    op.makeRHS(rhs, sourceTerms);
    return rhs;
  };

  FETimeDependent TDsolver(op.matrix(), rhs_f, grid);

  vector<double> internalNodes;

//   TDsolver.solve(10, .01, internalNodes, initial_conditions, rhs_f);
  TDsolver.solve_write(10, .01, internalNodes, initial_conditions, rhs_f, "vtk_output/solution");
  
  return 0;
  
}
