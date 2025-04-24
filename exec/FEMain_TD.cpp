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
  return 4.0*sin(5.0*Rsquared*sin(time));

}

float derivedf(double time, array<double, DIM> x)
{
  double Rsquared = (x[1]-9)*(x[1]-9)+x[0]*x[0];

  //calculations for partials:
  // dphi/dx = 4.0*cos(5.0*Rsquared*sin(time))*5.0*2.0*x[0]*sin(time)
  // d^2phi/dx^2 = -4.0*sin(5.0*Rsquared*sin(time))*5.0*2.0*x[0]*sin(time)*5.0*2.0*x[0]*sin(time)+4.0*cos(5.0*Rsquared*sin(time))*5.0*2.0*sin(time)
  
  // dphi/dy = 4.0*cos(5.0*Rsquared*sin(time))*5.0*2.0*(x[1]-9)*sin(time)
  // d^2phi/dy^2 = -4.0*sin(5.0*Rsquared*sin(time))*5.0*2.0*(x[1]-9)*sin(time)*5.0*2.0*(x[1]-9)*sin(time)+4.0*cos(5.0*Rsquared*sin(time))*5.0*2.0*sin(time)

  // dphi/dt = 4.0*cos(5.0*Rsquared*sin(time))*5.0*Rsquared*cos(time)
  double d2phidx2 = -4.0*sin(5.0*Rsquared*sin(time))*5.0*2.0*x[0]*sin(time)*5.0*2.0*x[0]*sin(time)+4.0*cos(5.0*Rsquared*sin(time))*5.0*2.0*sin(time);
  double d2phidy2 = -4.0*sin(5.0*Rsquared*sin(time))*5.0*2.0*(x[1]-9)*sin(time)*5.0*2.0*(x[1]-9)*sin(time)+4.0*cos(5.0*Rsquared*sin(time))*5.0*2.0*sin(time);
  double dphidt = 4.0*cos(5.0*Rsquared*sin(time))*5.0*Rsquared*cos(time);

  return d2phidx2 + d2phidy2 - dphidt;

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

  FEPoissonOperator<double> op(grid);
  

  vector<double> initial_conditions(grid.getNumNodes());

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

  function<double(const Node &, double)> boundary_cond = [](const Node &n, double time)
  {
    double x = n.getPosition()[0];
    double y = n.getPosition()[1];
    double Rsquared = (y-9)*(y-9)+x*x;
    if (Rsquared <= 26)
    {
      return -1.0;
    } else if (Rsquared >= 35)
    {
      return -1.0;
    }
    return 0.0;
  };

  FETimeDependent TDsolver(op.matrix(), rhs_f, grid);

  vector<double> phi_nodes;

  TDsolver.solve_write(1, .005, phi_nodes, initial_conditions, boundary_cond, "vtk_output/solution");



  //error analysis using infinity norm
  double maxerr = 0.0;

  for (int nodeidx = 0; nodeidx < grid.getNumNodes(); nodeidx ++)
  {
    Node n = grid.node(nodeidx);
    if (n.isInterior())
    {
      double newerr = abs(phi_nodes[n.getInteriorNodeID()] - sourcePhi(1, n.getPosition()));
      if (newerr > maxerr)
      {
        maxerr = newerr;
      }
    }
  }

  cout << "infinity norm of calculated phi and source phi is: " << maxerr << endl;
  
  return 0;
  
}
