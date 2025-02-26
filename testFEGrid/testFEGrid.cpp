#include <cstdio>   
#include <cmath>    
#include <cstring>  
#include <cassert>  
#include <vector>
#include "Node.H"   
#include "Element.H"
#include "FEGrid.H"

int main(int argc, char** argv)
{
  assert(argc == 2);
  std::string prefix(argv[1]);
  string nodeFile=prefix+".node";
  string eleFile =prefix+".ele";
  
  FEGrid grid(nodeFile, eleFile);
  
  vector<double> nodeData(grid.getNumNodes());
  array<double, DIM> p;
  for(unsigned int i=0; i<nodeData.size(); i++)
    {
      const Node& n = grid.node(i);
      p = n.getPosition();
      nodeData[i]=p[0];
    }
  FEWrite(&grid, &nodeData, "gridFile");

  return 0;
}

  
