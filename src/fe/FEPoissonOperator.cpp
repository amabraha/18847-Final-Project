#include <cmath>
#include <cassert>
#include <cmath>
#include <vector>
using namespace std;

#include "Node.H"
#include "Element.H"
#include "FEGrid.H"
#include "SparseMatrix.H"
#include "FEPoissonOperator.H"
#include <complex>

double dotprod(array<double, DIM> x, array<double, DIM> y)
{
    return x[0] * y[0] + x[1] * y[1];
}

template <typename T>
FEPoissonOperator<T>::FEPoissonOperator(const FEGrid &a_grid)
{
    m_grid = a_grid;
    // include all nodes in the matrix (both interior and boundary)
    m_matrix = SparseMatrix<T>(m_grid.getNumNodes(),  // rows
                               m_grid.getNumNodes()); // cols

    // loop through each element
    for (int iEl = 0; iEl < m_grid.getNumElts(); iEl++)
    {
        // access element with id iEl
        const Element &e = m_grid.element(iEl);

        // loop through the local vertices of element e
        for (int iVert = 0; iVert < VERTICES; ++iVert)
        {
            // 1) grab global node indices from the element connectivity
            int ig = e[iVert];                  // global row
            int jg = e[(iVert + 1) % VERTICES]; // global col

            // 2) diagonal entry L_{ig,ig} += area * grad_i · grad_i
            {
                array<int, 2> idx = {ig, ig};
                m_matrix[idx] += m_grid.elementArea(iEl) * dotprod(m_grid.gradient(iEl, iVert),
                                                                   m_grid.gradient(iEl, iVert));
            }

            // 3) off-diagonal entry L_{ig,jg} and symmetric partner L_{jg,ig}
            {
                array<int, 2> idx = {ig, jg};
                array<int, 2> idxT = {jg, ig};
                m_matrix[idx] += m_grid.elementArea(iEl) * dotprod(m_grid.gradient(iEl, iVert),
                                                                   m_grid.gradient(iEl, (iVert + 1) % VERTICES));
                m_matrix[idxT] = m_matrix[idx];
            }
        }
    }

    //apply boundary conditions to L
    for (int i = 0; i < m_grid.getNumNodes(); ++i)
    {
        const Node &n = m_grid.node(i);
        if (n.isInterior())
            continue;

        m_matrix.zeroRow(i); // wipe old Laplace row
        array<int, 2> d = {i, i};
        m_matrix[d] = T(1); // 1·φ_i = g
    }
}

template <typename T>
void FEPoissonOperator<T>::makeRHS(
    vector<T> &a_rhsAtNodes,
    const vector<T> &a_FCentroids) const
{
    // 1) global resize
    a_rhsAtNodes.resize(m_grid.getNumNodes());

    // 2) zero out
    for (int i = 0; i < (int)a_rhsAtNodes.size(); ++i)
        a_rhsAtNodes[i] = T(0.0);

    // 3) loop elements
    for (int iEl = 0; iEl < m_grid.getNumElts(); ++iEl)
    {
        T area = static_cast<T>(m_grid.elementArea(iEl));
        T factor = T(1.0) / T(3.0);

        const Element &e = m_grid.element(iEl);
        for (int iVert = 0; iVert < VERTICES; ++iVert)
        {
            int ig = e[iVert]; // global node index
            a_rhsAtNodes[ig] += area * a_FCentroids[iEl] * factor;
        }
    }
}

template <typename T>
void FEPoissonOperator<T>::applyDirichletBC(
    SparseMatrix<T> &L,
    std::vector<T> &rhs,
    std::function<T(const Node &)> Phi_omega) const
{
    // loop by *global* node index
    for (int i = 0; i < m_grid.getNumNodes(); ++i)
    {
        const Node &n = m_grid.node(i);
        if (n.isInterior())
            continue;

        T g = Phi_omega(n);

        L.zeroRow(i); // wipe old Laplace row
        array<int, 2> d = {i, i};
        L[d] = T(1); // 1·φ_i = g
        rhs[i] = g;  // enforce value
    }
}

template <typename T>
const FEGrid &FEPoissonOperator<T>::getFEGrid() const
{
    return m_grid;
}

template <typename T>
const SparseMatrix<T> &FEPoissonOperator<T>::matrix() const
{
    return m_matrix;
}

// Explicit template instantiations
template class FEPoissonOperator<float>;
template class FEPoissonOperator<double>;
template class FEPoissonOperator<std::complex<float>>;
template class FEPoissonOperator<std::complex<double>>;

