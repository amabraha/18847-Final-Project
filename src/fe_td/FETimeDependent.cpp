#include "FETimeDependent.H"
#include <array>
#include <cassert>
#include <iostream>
#include <cstring>
using namespace std;



FETimeDependent::FETimeDependent(const SparseMatrix& a_L, const function<vector<double>(double)>& a_f, const FEGrid& a_grid)
{
  m_L = a_L;
  m_f = a_f;
  m_grid = a_grid;
}

void FETimeDependent::step(double time, double dt, vector<double>& a_phi_out)
{
  vector<double> rhs = vector<double>(a_phi_out.size());
  SparseMatrix A(m_L, -1);

  for (int i = 0; i < a_phi_out.size(); i ++)
  {
    rhs[i] = a_phi_out[i]/dt - m_f(time+dt)[i];
    A.access(array<int, 2>{i,i}) += 1/dt;
  }

  JacobiSolver solver;
  //ehh just picked 1000 and 1e-6 arbitrarily
  solver.solve(a_phi_out, A, rhs, 1E-6, 1000);
}

void FETimeDependent::solve(double time, double dt, vector<double>& a_phi_out, 
                            const vector<double>& a_initial_cond, const function<vector<double>(double)>& a_boundary_cond)
{
  //initialize our output to the initial condition
  a_phi_out.resize(a_initial_cond.size());
  for (int i = 0; i < a_initial_cond.size(); i ++)
  {
    a_phi_out[i] = a_initial_cond[i];
  }

  double t = 0;
  //keep stepping by dt until we can't anymore
  while(t < time)
  {
    if (t + dt > time)
    {
      step(t, time-t, a_phi_out);
      break;
    }
    t += dt;
    step(t, dt, a_phi_out);
  }
}

void FETimeDependent::solve_write(double time, double dt, vector<double>& a_phi_out, 
                                  const vector<double>& a_initial_cond, const function<vector<double>(double)>& a_boundary_cond,
                                  string a_filename)
{
  //initialize our output to the initial condition
  a_phi_out.resize(a_initial_cond.size());
  for (int i = 0; i < a_initial_cond.size(); i ++)
  {
    a_phi_out[i] = a_initial_cond[i];
  }

  double t = 0;
  int file_counter = 0;
  //keep stepping by dt until we can't anymore
  while(t < time)
  {
    if (t + dt > time)
    {
      step(t, time-t, a_phi_out);
      t = time;
    } 
    else 
    {
      t += dt;
      step(t, dt, a_phi_out);
    }

    vector<double> phi_full;
    reinsert(m_grid, a_phi_out, phi_full);

    string filename = a_filename+to_string(file_counter)+".vtk";
    FEWrite(&m_grid, &phi_full, filename.c_str());
    file_counter += 1;

  }
}