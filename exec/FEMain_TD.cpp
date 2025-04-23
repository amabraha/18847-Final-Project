#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "FETimeDependent.H"
#include "ReInsert.H"
#include "JacobiSolver.H"
#include <iostream>
#include <array>
using namespace std;


float sourcePhi(double time, array<double, DIM> x)
{
  double Rsquared = (x[1]-9)*(x[1]-9)+x[0]*x[0];
  
  if (Rsquared <= 25) //inner boundary
  {
    return -1.0;
  } else if (Rsquared >= 36) //outer boundary
  {
    return 1.0;
  } else
  {
    return 4.0*sin(Rsquared*sin(time));
  }
}

float derivedf(double time, array<double, DIM> x)
{
  double Rsquared = (x[1]-9)*(x[1]-9)+x[0]*x[0];
  
  if (Rsquared <= 25) //inner boundary
  {
    return 0.0;
  } else if (Rsquared >= 36) //outer boundary
  {
    return 0.0;
  } else
  {
    // dphi/dx = 4.0*cos(Rsquared*sin(time))*2.0*x[0]*sin(time)
    // d^2phi/dx^2 = -4.0*sin(Rsquared*sin(time))*2.0*x[0]*sin(time)*2.0*x[0]*sin(time)+4.0*cos(Rsquared*sin(time))*2.0*sin(time)
    
    // dphi/dy = 4.0*cos(Rsquared*sin(time))*2.0*(x[1]-9)*sin(time)
    // d^2phi/dy^2 = -4.0*sin(Rsquared*sin(time))*2.0*(x[1]-9)*sin(time)*2.0*(x[1]-9)*sin(time)+4.0*cos(Rsquared*sin(time))*2.0*sin(time)

    // dphi/dt = 4.0*cos(Rsquared*sin(time))*Rsquared*cos(time)
    double d2phidx2 = -4.0*sin(Rsquared*sin(time))*2.0*x[0]*sin(time)*2.0*x[0]*sin(time)+4.0*cos(Rsquared*sin(time))*2.0*sin(time);
    double d2phidy2 = -4.0*sin(Rsquared*sin(time))*2.0*(x[1]-9)*sin(time)*2.0*(x[1]-9)*sin(time)+4.0*cos(Rsquared*sin(time))*2.0*sin(time);
    double dphidt = 4.0*cos(Rsquared*sin(time))*Rsquared*cos(time);

    return d2phidx2 + d2phidy2 - dphidt;
  }

}

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
  if(argc > 3)
      {
      cout << "this program takes one argument that is the .node and .ele OR .poly and max area";
      cout << "file prefix"<<endl;
      return 1;
    }
  FEGrid grid;
  string prefix(argv[1]);
  if (argc == 2)
  {
    string nodeFile = prefix+".node";
    string eleFile  = prefix+".ele";

    grid = FEGrid(nodeFile, eleFile);
  } else
  {
    string polyFile = prefix+".poly";
    double max_area = stod(argv[2]);
    
    grid = FEGrid(polyFile, max_area);
  }

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
    sourceTerms[i] = derivedf(time, centroid);
    }
  
    vector<double> rhs;
    op.makeRHS(rhs, sourceTerms);
    return rhs;
  };

  FETimeDependent TDsolver(op.matrix(), rhs_f, grid);

  vector<double> internalNodes;

//   TDsolver.solve(10, .01, internalNodes, initial_conditions, rhs_f);
  TDsolver.solve_write(4, .01, internalNodes, initial_conditions, rhs_f, "vtk_output/solution");
  
  return 0;
  
}
