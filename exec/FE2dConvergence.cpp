// exec/FE2dConvergence.cpp
// Convergence study for the 2‑D Poisson problem using P1 elements.
// Builds a sequence of Triangle meshes with successively smaller
// global area bound, solves –Δφ = f with exact Dirichlet data, and
// prints the observed order of accuracy.
//
// Usage  (compile with DIM=2):
//   fe2d_convergence.exe <poly‑prefix> <initial-max-area> [num-levels]
// -------------------------------------------------------------------
#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "JacobiSolver.H"
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

using std::array;
using std::cout;
using std::string;
using std::vector;

// analytic exact solution and forcing ----------------------------------------
static auto Phi_exact = [](const array<double, DIM> &X) -> double
{
    // smooth function that vanishes on the unit square boundary
    // return 2*X[0]+X[1];
    return X[0]*X[0]+X[1]*X[1];
};

static auto source2D = [](const array<double, DIM> &X) -> double
{
    // −Δ Φ = 2π² Φ for the choice above
    return -4.0;
};

// mesh‑size estimate ----------------------------------------------------------
static double estimate2Dh(const FEGrid &grid)
{
    double A = 0; // average triangle area
    for (int i = 0; i < grid.getNumElts(); i++) {
        A += grid.elementArea(i);
    }
    A /= grid.getNumElts();
    cout << "A: " << A << endl;
    return 2.0 * std::sqrt(A);      // edge length of iso‑area square
}

// RHS assembly helper ---------------------------------------------------------
static void build2DRHS(const FEGrid &grid,
                       const FEPoissonOperator<double> &op,
                       const decltype(source2D) &f,
                       vector<double> &rhs)
{
    int nel = grid.getNumElts();
    vector<double> fcent(nel);
    for (int e = 0; e < nel; ++e)
        fcent[e] = f(grid.centroid(e));

    op.makeRHS(rhs, fcent);
}

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
    if (argc < 3 || argc > 4)
    {
        cout << "Usage:\n  " << argv[0]
             << " <poly-prefix> <initial-max-area> [num-levels]\n";
        return 1;
    }
    string prefix = argv[1];        // ../FEData/test  (reads test.poly)
    double A0 = std::stod(argv[2]); // e.g. 0.005
    int levels = (argc == 4 ? std::stoi(argv[3]) : 5);

    vector<double> hs, errs;

    for (int L = 0; L < levels; ++L)
    {
        double max_area = A0 / std::pow(2, L);

        cout << "max_area: " << max_area << endl;

        // build 2‑D mesh via Triangle (no extrusion when DIM=2)
        FEGrid grid(prefix + ".poly", max_area);

        // assemble operator and RHS
        FEPoissonOperator<double> op(grid);
        vector<double> rhs;
        build2DRHS(grid, op, source2D, rhs);

        // copy matrix so we can modify it for Dirichlet BC
        SparseMatrix<double> A = op.matrix();
        op.applyDirichletBC(A, rhs, [](const Node &n)
                            { return Phi_exact(n.getPosition()); });

        // solve with Jacobi
        vector<double> phi;
        JacobiSolver<double> solver;
        solver.solve(phi, A, rhs, 1e-8, 100000);

        // compute nodal max error
        double maxErr = 0.0;
        for (int i = 0; i < grid.getNumNodes(); ++i)
        {
            maxErr = std::max(maxErr, std::abs(phi[i] - Phi_exact(grid.node(i).getPosition())));
        }

        double h = estimate2Dh(grid);
        hs.push_back(h);
        errs.push_back(maxErr);

        cout << "Level " << L
             << "   #tri=" << grid.getNumElts()
             << "   h=" << h
             << "   err=" << maxErr << "\n";
        string filename = "solution"+to_string(L)+".vtk";
        FEWrite(&grid, &phi, filename.c_str());
    }

    // print observed orders
    cout << "\nEstimated order p:\n";
    for (int i = 1; i < static_cast<int>(hs.size()); ++i)
    {
        double p = std::log(errs[i - 1] / errs[i]) / std::log(hs[i - 1] / hs[i]);
        cout << "  between lvl " << (i - 1) << "->" << i << " : p ≈ " << p << "\n";
    }
    return 0;
}
