#include "JacobiSolver.H"
#include <cassert>
#include <iostream>
#include <array>
using namespace std;

double norm(const vector<double>& a_v)
{
  double rtn = 0;
  for(unsigned int i=0; i<a_v.size(); i++)
    {
      double x = a_v[i];
      if(x < 0)
	{
	  if(rtn < -x) rtn=-x;
	}
      else
	{
	  if(rtn < x) rtn = x;
	}
    }
  return rtn;
}
