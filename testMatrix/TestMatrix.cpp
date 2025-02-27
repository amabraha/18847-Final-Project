#include "SparseMatrix.H"

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2)
    {
      cout << "this program takes one argument that is number of nodes\n ";
      return 1;
    }
  int cells = atoi(argv[1]);

  if(cells< 4)
    {
      cout<<"you'll need at least 4 nodes to apply the Laplacian \n";
      return 2;
    }

  SparseMatrix A(cells-2, cells);
  array<int, 2> left   = {{0,0}};
  array<int, 2> center = {{0,1}};
  array<int, 2> right  = {{0,2}};
  for(int i=0; i<cells-2; i++, left[0]++, left[1]++, center[0]++, center[1]++,
	right[0]++, right[1]++)
    {
      A[left]=1;
      A[center]=-2;
      A[right]=1;
    }

  //  now that A is built, let's try some simple functions
  // first, linear function. The second derivative is zero, and doubling point can 
  // represent these small integer values and their multiples by 1 and 2 exactly.

  vector<double> phi(cells);
  for(int i=0; i<cells; i++)
    {
      phi[i] = i;
    }

  vector<double> LOfPhi = A*phi;
  bool pass = true;
  for(int i=0; i<cells-2; i++)
    {
      if(LOfPhi[i] != 0)
	{
	  cout<< "error in matrix multiply of linear phi data LOfPhi[";
	  cout<<i<<"] = "<<LOfPhi[i]<<endl;
	  pass=false;
	}
    }
  if(pass)
    {
      cout<<"linear phi passed\n";
    }
  else
    {
      cout<<"linear phi failed\n";
    }
  
  // now a quadratic function, also exactly representable in this range of doubling point
  int j=-cells/2;
  for(int i=0; i<cells; j++, i++)
    {
      phi[i] = j*j - i;
    }

  LOfPhi = A*phi;

  pass = true;
  for(int i=0; i<cells-2; i++)
    {
      if(LOfPhi[i] != 2)
	{
	  cout<< "error in matrix multiply of linear phi data LOfPhi[";
	  cout<<i<<"] = "<<LOfPhi[i]<<endl;
	  pass=false;
	}
    }
  if(pass)
    {
      cout<<"quadratic phi passed\n";
    }
  else
    {
      cout<<"quadratic phi failed\n";
    }
  
  cout<<endl;

  
  return 0;
}
