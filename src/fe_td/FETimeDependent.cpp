#include "FETimeDependent.H"
#include <array>
#include <cassert>
#include <iostream>
using namespace std;



FETimeDependent::FETimeDependent(const SparseMatrix& a_L, const function<vector<double>(double)>& a_f)
{
    m_L = a_L;
    m_f = a_f;
}

void FETimeDependent::step(double time, double dt, vector<double>& a_phi_out)
{
    vector<double> rhs = vector<double>(a_phi_out.size());
    SparseMatrix A(m_L);

    for (int i = 0; i < a_phi_out.size(); i ++)
    {
        rhs[i] = - a_phi_out[i]/dt + m_f(time+dt)[i];
        A.access(array<int, 2>{i,i}) -= 1/dt;
    }

    JacobiSolver solver;
    solver.solve(a_phi_out, A, rhs, 100, 100);
}

void FETimeDependent::solve(double time, double dt, vector<double>& a_phi_out, 
                       const vector<double>& a_initial_cond, const function<vector<double>(double)>& a_boundary_cond)
{
    //initialize our output to the initial condition
    for (int i = 0; i < a_initial_cond.size(); i ++)
    {
        a_phi_out[i] = a_initial_cond[i];
    }

    double t = 0;
    while(t <= time)
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