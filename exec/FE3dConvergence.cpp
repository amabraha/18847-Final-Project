// exec/FE3dConvergence.cpp
#include "FEGrid.H"
#include "FEPoissonOperator.H"
#include "JacobiSolver.H"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// -----------------------------------------------------------------------------
// Exact analytic solution Φ  and forcing  f = -ΔΦ
static auto Phi_exact = [](const array<double, DIM> &X) -> double
{
    return sin(M_PI * X[0]) * sin(M_PI * X[1]) * sin(M_PI * X[2]);
};
static auto source3D = [](const array<double, DIM> &X) -> double
{
    return 3.0 * M_PI * M_PI * Phi_exact(X); // -ΔΦ
};

// -----------------------------------------------------------------------------
// Mesh spacing  h = 2 * (avg tetra volume)^{1/3}
static double estimate3Dh(const FEGrid &grid)
{
    double Vtot = 0.0;
    for (int e = 0; e < grid.getNumElts(); ++e)
        Vtot += grid.elementArea(e);
    double Vavg = Vtot / grid.getNumElts();
    return 2.0 * cbrt(Vavg);
}

// -----------------------------------------------------------------------------
// Build global RHS from  f  via centroid sampling
static void build3DRHS(const FEGrid &grid,
                       FEPoissonOperator<double> &op,
                       const decltype(source3D) &f,
                       vector<double> &rhs)
{
    int nEl = grid.getNumElts();
    vector<double> fcent(nEl);
    for (int i = 0; i < nEl; ++i)
        fcent[i] = f(grid.centroid(i));
    op.makeRHS(rhs, fcent);
}

// -----------------------------------------------------------------------------
// Main driver
int main(int argc, char **argv)
{
    if (argc < 3 || argc > 4)
    {
        cout << "Usage:\n  " << argv[0]
             << " <poly_prefix> <initial_max_area> [num_levels]\n";
        return 1;
    }

    string prefix = argv[1];   // ../FEData/test
    double A0 = stod(argv[2]); // e.g. 0.01
    int levels = (argc == 4 ? stoi(argv[3]) : 4);

    vector<double> hs, errors;

    for (int L = 0; L < levels; ++L)
    {
        double max_area = A0 / pow(2.0, L);

        FEGrid grid(prefix + ".poly", max_area); // 2D→3D extrusion

        FEPoissonOperator<double> op(grid);
        SparseMatrix<double> A = op.matrix();

        vector<double> rhs;
        build3DRHS(grid, op, source3D, rhs);

        op.applyDirichletBC(
            A, rhs,
            [&](const Node &n)
            { return Phi_exact(n.getPosition()); });

        vector<double> phi;
        JacobiSolver<double> solver;
        solver.solve(phi, A, rhs, 1e-8, 5000);

        // infinity-norm error
        double maxErr = 0.0;
        for (int i = 0; i < grid.getNumNodes(); ++i)
            maxErr = max(maxErr, fabs(phi[i] - Phi_exact(grid.node(i).getPosition())));

        double h = estimate3Dh(grid);
        hs.push_back(h);
        errors.push_back(maxErr);

        cout << "Level " << L
             << "   #tets=" << grid.getNumElts()
             << "   h=" << h
             << "   err=" << maxErr << '\n';
    }

    // convergence orders
    cout << "\nEstimated order p:\n";
    const double eps = 1e-12;
    for (int i = 1; i < (int)hs.size(); ++i)
    {
        double ratio = hs[i - 1] / hs[i];
        if (fabs(ratio - 1.0) < eps)
            cout << "  between lvl " << i - 1 << "→" << i
                 << " : h unchanged (refine further)\n";
        else
        {
            double p = log(errors[i - 1] / errors[i]) / log(ratio);
            cout << "  between lvl " << i - 1 << "→" << i
                 << " : p ≈ " << p << '\n';
        }
    }
    return 0;
}
