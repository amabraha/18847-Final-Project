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

const SparseMatrix& FEPoissonOperator::matrix() const
{
  return m_matrix;
}
