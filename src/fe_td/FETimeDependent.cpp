#include "FETimeDependent.H"
#include <array>
#include <cassert>
#include <iostream>
#include <cstring>
using namespace std;



FETimeDependent::FETimeDependent(const SparseMatrix<double>& a_L, const function<vector<double>(double)>& a_f, const FEGrid& a_grid)
{
  m_L = a_L;
  m_f = a_f;
  m_grid = a_grid;
}

void FETimeDependent::apply_boundary(vector<double>& a_phi, double time, const function<double(const Node &, double)>& a_boundary_cond)
{
  for (int i = 0; i < m_grid.getNumNodes(); ++i)
  {
      const Node &n = m_grid.node(i);
      if (n.isInterior())
          continue;

      double g = a_boundary_cond(n, time);
      a_phi[i] = g;  // enforce value
  }
}

void FETimeDependent::step(double time, double dt, vector<double>& a_phi_out, const function<double(const Node &, double)>& a_boundary_cond)
{
  apply_boundary(a_phi_out, time, a_boundary_cond);

  vector<double> rhs = vector<double>(a_phi_out.size());
  //copy -m_L to A
  SparseMatrix<double> A(m_L, 1);

  for (int i = 0; i < a_phi_out.size(); i ++)
  {
    //for the rhs this constructs 1/dt phi(t) - f(t + dt)
    rhs[i] = a_phi_out[i]/dt + m_f(time+dt)[i];

    //this constructs 1/dt I - L
    A.access(array<int, 2>{i,i}) += 1/dt;
  }

  JacobiSolver<double> solver;
  //ehh just picked 1000 and 1e-6 arbitrarily
  solver.solve(a_phi_out, A, rhs, 1E-7, 1000);

  apply_boundary(a_phi_out, time+dt, a_boundary_cond);
}

void FETimeDependent::solve(double time, double dt, vector<double>& a_phi_out, 
                            const vector<double>& a_initial_cond, const function<double(const Node &, double)>& a_boundary_cond)
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
      step(t, time-t, a_phi_out, a_boundary_cond);
      break;
    }
    t += dt;
    step(t, dt, a_phi_out, a_boundary_cond);
  }
}

void FETimeDependent::solve_write(double time, double dt, vector<double>& a_phi_out, 
                                  const vector<double>& a_initial_cond, const function<double(const Node &, double)>& a_boundary_cond,
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
      step(t, time-t, a_phi_out, a_boundary_cond);
      t = time;
    } 
    else 
    {
      t += dt;
      step(t, dt, a_phi_out, a_boundary_cond);
    }

    //visit writing process

    string filename = a_filename+to_string(file_counter)+".vtk";
    FEWrite(&m_grid, &a_phi_out, filename.c_str());
    file_counter += 1;

  }
}