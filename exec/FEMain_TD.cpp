#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "FETimeDependent.H"
#include "ReInsert.H"
#include "JacobiSolver.H"
#include <iostream>
#include <array>
using namespace std;


double sourcePhi(double time, array<double, DIM> x)
{
  return 2*x[0]+x[1]+ 3*time;

}

double derivedf(double time, array<double, DIM> x)
{
  //d2phidx2 + d2phidy2 + dphidt;
  return 0+3;
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

  for (int nodeidx = 0; nodeidx < grid.getNumNodes(); nodeidx ++)
    {
      Node n = grid.node(nodeidx);
      initial_conditions[nodeidx] = sourcePhi(0, n.getPosition());
    }

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
    return sourcePhi(time, n.getPosition());
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

  double finaltime = 1;
  double timestep = .005;

  TDsolver.solve_write(finaltime, timestep, phi_nodes, initial_conditions, boundary_cond, "vtk_output/solution");



  //error analysis using infinity norm
  double maxerr = 0.0;

  for (int nodeidx = 0; nodeidx < grid.getNumNodes(); nodeidx ++)
  {
    Node n = grid.node(nodeidx);
    if (n.isInterior())
    {
      double newerr = abs(phi_nodes[nodeidx] - sourcePhi(finaltime, n.getPosition()));
      if (newerr > maxerr)
      {
        maxerr = newerr;
      }
    }
  }
  cout << "infinity norm of difference between calculated phi and source phi is: " << maxerr << endl;
 
   
 
  double t = 0;
  double dt = timestep;
  int file_counter = 0;


  vector<double> phi_vector(grid.getNumNodes());
  //keep stepping by dt until we can't anymore
  while(t < finaltime)
  {
    if (t + dt > finaltime)
    {
      t = finaltime;
    } 
    else 
    {
      t += dt;
    }
    
  

    for (int nodeidx = 0; nodeidx < grid.getNumNodes(); nodeidx ++)
    {
      Node n = grid.node(nodeidx);
      phi_vector[nodeidx] = sourcePhi(t, n.getPosition());
    }

    //visit writing process

    string filename = string("phi_output/solution")+to_string(file_counter)+".vtk";
    FEWrite(&grid, &phi_vector, filename.c_str());
    file_counter += 1;

  }
  
  return 0;
  
}
