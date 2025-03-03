#include "SparseMatrix.H"
#include "JacobiSolver.H"

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2)
    {
      cout << "this program takes one argument that is the dimension of the matrix\n ";
      return 1;
    }
  int N = atoi(argv[1]);

  if(N < 1)
    {
      cout<<"Dimension of matrix should be at least 1\n";
      return 2;
    }

  // we are going to construct a tridiagonal matrix with diagonal entries
  // 3 on the diagonal. 1 on the upper off diagonal. 2 on the lower off diagonal.
  SparseMatrix A(N, N);
  

  array<int, 2> diagonal = array<int, 2>{0,0};
  A[diagonal] = 3.0;

  diagonal = array<int, 2>{1,1};
  array<int, 2> upper_diagonal = array<int, 2>{0,1};
  array<int, 2> lower_diagonal = array<int, 2>{1,0};

  for(int i=0; i<N-1; i++, diagonal[0]++, diagonal[1]++, upper_diagonal[0]++, upper_diagonal[1]++,
	lower_diagonal[0]++, lower_diagonal[1]++)
    {
      A[diagonal]=3.0;
      A[upper_diagonal]=1.0;
      A[lower_diagonal]=2.0;
    }

  // make random vector v. We will compute the solution to Av=b
  vector<double> v = vector<double>(N);
  for(int i=0; i < N; i++)
    {
      v[i] = i;
    }
  vector<double> b = A*v;

  //vref will be the solution for v given by our Jacobi solver
  vector<double> vref = vector<double>(N);
  JacobiSolver JS = JacobiSolver();
  JS.solve(vref, A, b, .1, N);
  bool pass = true;
  for(int i=0; i<N; i++)
    {
      if(abs(vref[i]-v[i]) > .2 && abs(vref[i]-v[i]) > .2*v[i])
	{
    printf("Jacobi iteration solution at %d is too off\nexpected v[%d] = %f, got back vref[%d] = %f\n", i, i, v[i], i, vref[i]);
	  pass=false;
	}
    }
  if(pass)
    {
      cout<<"Jacobi iteration passed\n";
    }
  else
    {
      cout<<"Jacobi iteration failed\n";
    }
  cout << endl;

  
  return 0;
}
