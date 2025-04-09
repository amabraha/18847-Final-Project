#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <vector>
#include <cassert>
#include <cmath>
#include <array>
using namespace std;
class SparseMatrix
{
public:
  /// set up an M rows and N columns sparse matrix with all values of zero (no non-zero elements)
  SparseMatrix();
  SparseMatrix(int a_M, int a_N);

  /// copy constructor
  SparseMatrix(const SparseMatrix &M);

  /// copy constructor with a multiplier of k
  SparseMatrix(const SparseMatrix &M, double k);

  /// Matrix Vector multiply.  a_v.size()==N, returns vector of size M
  /** 
      Important: you should be able to perform matrix-vector multiplication by
      iterating through the vectors defining the rows of the matrix.
      <code> (A*v)[i] </code> = \f$ \sum_j \f$ <code> m\_data[i][j] v[colIndex[i][j]] </code>.
      Note that the ordering of the!! columm elements of row i implied by colIndex[i]
      Does not have to be in any particular order.
   */
  vector<double> operator*(const vector<double>& a_v) const;

  ///accessor functions for get and set operations of matrix elements. 
  /**
    Non-const version of (*this)[p] returns a reference to the (p[0],p[1])
    element of the matrix. If there is no element of the matrix corresponding to
    that pair of indices, create a new element of the array before returning the 
    reference to the value. Note that you will have to search through the p[0] 
    row of the array to find the element. You should be able to create a new 
    element of the array using push_back, in light of the comment above.
  */
  double& operator[](array<int, 2>& a_index);
  double& access(array<int, 2> a_index);

  ///accessor function just to get a value
  const double& operator[](array<int, 2>& a_index) const;
  const double& access(array<int, 2> a_index) const;

  /// zero out all the elements, but leave the sparse structure in place.
  void zero();

  /// build and return a new SparseMatrix that is the transpose of the input matrix.
  SparseMatrix transpose() const;
  /// Returns the value of m_m.
  unsigned int M() const;
  /// Returns the value of m_n.
  unsigned int N() const;
  /// Checks to see whether the matrix is symmetric.
  bool symmetric() const;

  void print() const;
private:
  unsigned int m_m, m_n;
  double m_zero;
  /* 
  m_data[p][q]) contains the 
  (p,colIndex[p][q]) element of the array.
  */ 
  vector<vector<double> > m_data;
  vector<vector<int> >   m_colIndex;
};

#endif
