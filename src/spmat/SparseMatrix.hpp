#ifndef _SPARSEMATRIX_HPP_
#define _SPARSEMATRIX_HPP_

#include "SparseMatrix.H"
#include <complex>

template <typename T>
SparseMatrix<T>::SparseMatrix()
    : m_m(0), m_n(0), m_zero(T(0))
{
}

template <typename T>
SparseMatrix<T>::SparseMatrix(int a_M, int a_N)
{
    m_m = a_M;
    m_n = a_N;
    m_zero = T(0);
    m_data = vector<vector<T>>(a_M, vector<T>(0));
    m_colIndex = vector<vector<int>>(a_M, vector<int>(0));
}

template <typename T>
vector<T> SparseMatrix<T>::operator*(const vector<T> &a_v) const
{
    vector<T> res = vector<T>(m_m, T(0));

    for (int row = 0; row < m_m; row++)
    {
        for (int col = 0; col < m_colIndex[row].size(); col++)
        {
            res[row] += m_data[row][col] * a_v[m_colIndex[row][col]];
        }
    }

    return res;
}

template <typename T>
T &SparseMatrix<T>::operator[](array<int, 2> &a_index)
{
    int row = a_index[0];
    int col = a_index[1];
    int col_index = 0;

    for (; col_index < m_colIndex[row].size(); col_index++)
    {
        if (m_colIndex[row][col_index] == col)
        {
            return m_data[row][col_index];
        }
    }

    m_colIndex[row].push_back(col);
    m_data[row].push_back(T(0));

    return m_data[row][col_index];
}

template <typename T>
const T &SparseMatrix<T>::operator[](array<int, 2> &a_index) const
{
    int row = a_index[0];
    int col = a_index[1];
    int col_index = 0;

    for (; col_index < m_colIndex[row].size(); col_index++)
    {
        if (m_colIndex[row][col_index] == col)
        {
            return m_data[row][col_index];
        }
    }

    return m_zero;
}

template <typename T>
void SparseMatrix<T>::zero()
{
    for (int row = 0; row < m_m; row++)
    {
        for (int col = 0; col < m_colIndex[row].size(); col++)
        {
            m_data[row][col] = T(0);
        }
    }
}

template <typename T>
int SparseMatrix<T>::M() const
{
    return m_m;
}

template <typename T>
int SparseMatrix<T>::N() const
{
    return m_n;
}

template <typename T>
bool SparseMatrix<T>::symmetric() const
{
    for (int i = 0; i < m_m; i++)
    {
        for (int j = 0; j < m_colIndex[i].size(); j++)
        {
            int col = m_colIndex[i][j];
            array<int, 2> index1 = {i, col};
            array<int, 2> index2 = {col, i};
            if (m_data[i][j] != (*this)[index2])
            {
                return false;
            }
        }
    }
    return true;
}

// Explicit template instantiations
template class SparseMatrix<float>;
template class SparseMatrix<double>;
template class SparseMatrix<complex<float>>;
template class SparseMatrix<complex<double>>;

#endif