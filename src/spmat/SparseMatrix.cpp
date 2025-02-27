#include "SparseMatrix.H"
#include <iostream>

SparseMatrix::SparseMatrix()
  :m_m(0), m_n(0), m_zero(0)
{
}

SparseMatrix::SparseMatrix(int a_M, int a_N)
{
  m_m = a_M;
  m_n = a_N;
  m_zero = 0.0;
  m_data = vector<vector<double>>(a_M, vector<double>(0));
  m_colIndex = vector<vector<int>>(a_M, vector<int>(0));
}

vector<double> SparseMatrix::operator*(const vector<double>& a_v) const
{
  vector<double> res = vector<double>(m_m, 0.0);

  for (int row = 0; row < m_m; row ++)
    {
      for (int col = 0; col < m_colIndex[row].size(); col ++)
        {
          res[row] += m_data[row][col] * a_v[m_colIndex[row][col]];
        }
    }
  
    return res;
}

double& SparseMatrix::operator[](array<int, 2>& a_index)
{
  int row = a_index[0];
  int col = a_index[1];
  int col_index = 0;

  for (; col_index < m_colIndex[row].size(); col_index ++)
    {
      //search if row,col exists as a nonzero entry in our sparse matrix
      if (m_colIndex[row][col_index] == col)
        {
          return m_data[row][col_index];
        }
    }

  //if we can't find it, make a new one
  m_colIndex[row].push_back(col);
  m_data[row].push_back(0.0);

  return m_data[row][col_index];
}

const double& SparseMatrix::operator[](array<int, 2>& a_index) const
{
  int row = a_index[0];
  int col = a_index[1];
  int col_index = 0;

  for (; col_index < m_colIndex[row].size(); col_index ++)
    {
      //search if row,col exists as a nonzero entry in our sparse matrix
      if (m_colIndex[row][col_index] == col)
        {
          return m_data[row][col_index];
        }
    }

  //you shouldn't be trying to access a const value of an entry that doesn't exist.
  return m_data[row][col_index];
}

