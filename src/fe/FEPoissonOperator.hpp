#ifndef _FEPOISSONOPERATOR_HPP_
#define _FEPOISSONOPERATOR_HPP_

#include "FEPoissonOperator.H"

template <typename T>
T dotprod(array<T, DIM> x, array<T, DIM> y)
{
    return x[0] * y[0] + x[1] * y[1];
}

template <typename T>
FEPoissonOperator<T>::FEPoissonOperator(const FEGrid &a_grid)
{
    m_grid = a_grid;
    m_matrix = SparseMatrix<T>(m_grid.getNumInteriorNodes(), m_grid.getNumInteriorNodes());

    // loop through each element
    for (int iEl = 0; iEl < m_grid.getNumElts(); iEl++)
    {
        // access element with id iEl
        const Element &e = m_grid.element(iEl);

        // loop through the local vertices of element e
        for (int iVert = 0; iVert < VERTICES; iVert++)
        {
            // xn and xm are a pair of adjacent vertices in element e
            const Node &xn = m_grid.node(e[iVert]);
            const Node &xm = m_grid.node(e[(iVert + 1) % VERTICES]);

            // we only populate L_{n,n} if xn is interior
            if (xn.isInterior())
            {
                array<int, 2> sparse_matrix_indx = array<int, 2>{xn.getInteriorNodeID(), xn.getInteriorNodeID()};
                m_matrix[sparse_matrix_indx] += m_grid.elementArea(iEl) * dotprod<T>(m_grid.gradient(iEl, iVert), m_grid.gradient(iEl, iVert));
            }

            // we only populate L_{n,m} if xn and xm are interior
            if (xn.isInterior() && xm.isInterior())
            {
                array<int, 2> sparse_matrix_indx = array<int, 2>{xn.getInteriorNodeID(), xm.getInteriorNodeID()};
                array<int, 2> sparse_matrix_indx_trans = array<int, 2>{sparse_matrix_indx[1], sparse_matrix_indx[0]};

                m_matrix[sparse_matrix_indx] += m_grid.elementArea(iEl) * dotprod<T>(m_grid.gradient(iEl, iVert), m_grid.gradient(iEl, (iVert + 1) % VERTICES));
                m_matrix[sparse_matrix_indx_trans] = m_matrix[sparse_matrix_indx];
            }
        }
    }
}

template <typename T>
void FEPoissonOperator<T>::makeRHS(
    vector<T> &a_rhsAtNodes,
    const vector<T> &a_FCentroids) const
{
    a_rhsAtNodes.resize(m_grid.getNumInteriorNodes());

    // zero out all entries in RHS
    for (int i = 0; i < a_rhsAtNodes.size(); i++)
    {
        a_rhsAtNodes[i] = T(0);
    }

    // loop through each element
    for (int iEl = 0; iEl < m_grid.getNumElts(); iEl++)
    {
        // access element with id iEL
        const Element &e = m_grid.element(iEl);

        // loop through the local vertices of element e
        for (int iVert = 0; iVert < VERTICES; iVert++)
        {
            // xn is the corresponding vertex in element e
            const Node &xn = m_grid.node(e[iVert]);

            // we only populate b_n if xn is interior
            if (xn.isInterior())
            {
                a_rhsAtNodes[xn.getInteriorNodeID()] += m_grid.elementArea(iEl) * a_FCentroids[iEl] * T(1.0 / 3.0);
            }
        }
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

#endif