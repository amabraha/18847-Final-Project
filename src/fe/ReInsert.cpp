
#include "ReInsert.H"
#include "FEGrid.H"
#include <vector>

void reinsert(const FEGrid& a_grid, const vector<double>& a_internalSolution, vector<double>& a_fullVector)
{
  int nodes = a_grid.getNumNodes();
  a_fullVector.resize(0);
  a_fullVector.resize(nodes,0.0);

  for(int i=0; i<nodes; i++)
    {
      const Node& n = a_grid.node(i);
      if(n.isInterior())
	{
	  a_fullVector[i] = a_internalSolution[n.getInteriorNodeID()];
	}
    }
}
