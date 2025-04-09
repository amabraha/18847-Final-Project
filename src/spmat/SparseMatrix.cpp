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

SparseMatrix::SparseMatrix(const SparseMatrix &M)
{
  m_m = M.m_m;
  m_n = M.m_n;
  m_zero = 0.0;
  m_data = M.m_data;
  m_colIndex = M.m_colIndex;
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

double& SparseMatrix::access(array<int, 2> a_index)
{
  return (*this)[a_index];
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

  //since we aren't modifying our sparse matrix, just return 0.0
  return m_zero;
}

const double& SparseMatrix::access(array<int, 2> a_index) const
{
  return (*this)[a_index];
}

void SparseMatrix::zero()
{

  for (int row = 0; row < m_m; row ++)
    {
      for (int col = 0; col < m_colIndex[row].size(); col ++)
        {
          m_data[row][col] = 0.0;
        }
    }

}

SparseMatrix SparseMatrix::transpose() const
{

  SparseMatrix transposed = SparseMatrix(m_n, m_m);

  for (int row = 0; row < m_m; row ++)
    {
      for (int col = 0; col < m_colIndex[row].size(); col ++)
        {
          array<int, 2> transposed_index = array<int, 2>{m_colIndex[row][col], row};
          transposed[transposed_index] = m_data[row][col];
        }
    }

  return transposed;
}

unsigned int SparseMatrix::M() const
{
  return m_m;
}

unsigned int SparseMatrix::N() const
{
  return m_m;
}

bool SparseMatrix::symmetric() const
{

  for (int row = 0; row < m_m; row ++)
    {
      for (int col = 0; col < m_colIndex[row].size(); col ++)
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



void SparseMatrix::print() const
{
  printf("[");

  for (int row = 0; row < m_m; row ++)
    {
      for (int col = 0; col < m_colIndex[row].size(); col ++)
        {
          printf("(%d, %d, %f), ", row, m_colIndex[row][col], m_data[row][col]);
        }
    }
  printf("]\n");
}


