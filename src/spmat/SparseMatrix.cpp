#include "SparseMatrix.H"
#include <cassert>
#include <iostream>
#include <vector>
#include <complex>
#include <type_traits>
#include <cstdio>
using namespace std;

template <typename T>
SparseMatrix<T>::SparseMatrix()
    : m_m(0), m_n(0), m_zero(T(0))
{
}

template <typename T>
SparseMatrix<T>::SparseMatrix(int a_m, int a_n)
    : m_m(a_m), m_n(a_n), m_zero(T(0))
{
  m_data.resize(a_m);
  m_colIndex.resize(a_m);
}


template <typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T> &M)
{
  m_m = M.m_m;
  m_n = M.m_n;
  m_zero = T(0);
  m_data = M.m_data;
  m_colIndex = M.m_colIndex;
}


template <typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T> &M, T k)
{
  m_m = M.m_m;
  m_n = M.m_n;
  m_zero = T(0);
  m_colIndex = M.m_colIndex;
  m_data = M.m_data;

  for (int row = 0; row < m_m; row ++)
    {
      for (int col = 0; col < m_colIndex[row].size(); col ++)
        {
          m_data[row][col] = M.m_data[row][col] * k;
        }
    }
}

template <typename T>
vector<T> SparseMatrix<T>::operator*(const vector<T> &a_v) const
{
  vector<T> res = vector<T>(m_m, T(0.0));

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
    // search if row,col exists as a nonzero entry in our sparse matrix
    if (m_colIndex[row][col_index] == col)
    {
      return m_data[row][col_index];
    }
  }

  // if we can't find it, make a new one
  m_colIndex[row].push_back(col);
  m_data[row].push_back(T(0.0));

  return m_data[row][col_index];
}

template <typename T>
T& SparseMatrix<T>::access(array<int, 2> a_index)
{
  return (*this)[a_index];
}

template <typename T>
const T &SparseMatrix<T>::operator[](array<int, 2> &a_index) const
{
  int row = a_index[0];
  int col = a_index[1];
  int col_index = 0;

  for (; col_index < m_colIndex[row].size(); col_index++)
  {
    // search if row,col exists as a nonzero entry in our sparse matrix
    if (m_colIndex[row][col_index] == col)
    {
      return m_data[row][col_index];
    }
  }

  // since we aren't modifying our sparse matrix, just return 0.0
  return m_zero;
}

template <typename T>
const T& SparseMatrix<T>::access(array<int, 2> a_index) const
{
  return (*this)[a_index];
}

template <typename T>
void SparseMatrix<T>::zero()
{
  for (int row = 0; row < m_m; row++)
  {
    for (int col = 0; col < m_colIndex[row].size(); col++)
    {
      m_data[row][col] = T(0.0);
    }
  }
}

template <typename T>
SparseMatrix<T> SparseMatrix<T>::transpose() const
{

  SparseMatrix<T> transposed = SparseMatrix<T>(m_n, m_m);

  for (int row = 0; row < m_m; row++)
  {
    for (int col = 0; col < m_colIndex[row].size(); col++)
    {
      array<int, 2> transposed_index = array<int, 2>{m_colIndex[row][col], row};
      transposed[transposed_index] = m_data[row][col];
    }
  }

  return transposed;
}

template <typename T>
unsigned int SparseMatrix<T>::M() const
{
  return m_m;
}

template <typename T>
unsigned int SparseMatrix<T>::N() const
{
  return m_n;
}

template <typename T>
bool SparseMatrix<T>::symmetric() const
{
  for (int row = 0; row < m_m; row++)
  {
    for (int col = 0; col < m_colIndex[row].size(); col++)
    {
      array<int, 2> index = array<int, 2>{row, col};
      array<int, 2> transposed_index = array<int, 2>{col, row};
      if ((*this)[index] != (*this)[transposed_index])
      {
        return false;
      }
    }
  }

  return true;
}

// Overload #1: real types
template <typename U>
typename std::enable_if<std::is_floating_point<U>::value>::type
printValue(const U &x)
{
  printf("%f", static_cast<double>(x));
}

// Overload #2: complex types
template <typename U>
void printValue(const std::complex<U> &z)
{
  printf("%f + %fi",
         static_cast<double>(z.real()),
         static_cast<double>(z.imag()));
}
// --------------------------------------------------------

template <typename T>
void SparseMatrix<T>::print() const
{
  printf("[");
  for (int row = 0; row < m_m; ++row)
  {
    for (size_t ci = 0; ci < m_colIndex[row].size(); ++ci)
    {
      int col = m_colIndex[row][ci];
      printf("(%d, %d, ", row, col);

      // call the right overload at compile time
      printValue(m_data[row][ci]);

      printf("), ");
    }
  }
  printf("]\n");
}

/* === new helper body ============================================ */
template <typename T>
void SparseMatrix<T>::zeroRow(int r)
{
  assert(r >= 0 && r < static_cast<int>(m_m));
  m_colIndex[r] = vector<int>(0);
  m_data[r] = vector<T>(0);
}

// Explicit template instantiations
template class SparseMatrix<float>;
template class SparseMatrix<double>;
template class SparseMatrix<std::complex<float>>;
template class SparseMatrix<std::complex<double>>;
