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


double dotprod(array<double, DIM> x, array<double, DIM> y)
{
  return x[0]*y[0] + x[1]*y[1];
}

FEPoissonOperator::FEPoissonOperator(const FEGrid& a_grid)
{
  m_grid = a_grid;
  m_matrix = SparseMatrix(m_grid.getNumInteriorNodes(), m_grid.getNumInteriorNodes());

  //loop through each element
  for(int iEl = 0; iEl < m_grid.getNumElts(); iEl ++)
    {
      //access element with id iEl
      const Element& e = m_grid.element(iEl);

      //loop through the local vertices of element e
      for(int iVert = 0; iVert < VERTICES; iVert ++)
        {
          //xn and xm are a pair of adjacent vertices in element e
          const Node& xn = m_grid.node(e[iVert]);
          const Node& xm = m_grid.node(e[(iVert+1) % VERTICES]);

          //we will separately populate matrix entries of the form L_{n,n} and entries of the form L_{m,n} (where m,n are distinct)

          //we only populate L_{n,n} if xn is interior
          if (xn.isInterior())
            {
              array<int, 2> sparse_matrix_indx = array<int, 2>{xn.getInteriorNodeID(), xn.getInteriorNodeID()};

              //using the pseudocode in page 21 of lecture 8, we know the gradients are constant,
              // so integrating over Ke is the same as multiplying with the area of Ke
              m_matrix[sparse_matrix_indx] += m_grid.elementArea(iEl)*dotprod(m_grid.gradient(iEl, iVert), m_grid.gradient(iEl, iVert));
            }

          //we only populate L_{n,m} if xn and xm are interior
          if (xn.isInterior() && xm.isInterior())
            {
              //we will simultaneously update L_{n,m} and L_{m,n} since L is symmetric
              array<int, 2> sparse_matrix_indx = array<int, 2>{xn.getInteriorNodeID(), xm.getInteriorNodeID()};
              array<int, 2> sparse_matrix_indx_trans = array<int, 2>{sparse_matrix_indx[1], sparse_matrix_indx[0]};

              //using the pseudocode in page 21 of lecture 8, we know the gradients are constant,
              // so integrating over Ke is the same as multiplying with the area of Ke
              m_matrix[sparse_matrix_indx] += m_grid.elementArea(iEl)*dotprod(m_grid.gradient(iEl, iVert), m_grid.gradient(iEl, (iVert+1) % VERTICES));
              m_matrix[sparse_matrix_indx_trans] = m_matrix[sparse_matrix_indx];
            }
        }
    }
}

void FEPoissonOperator::makeRHS(
  vector<double> & a_rhsAtNodes, 
  const vector<double> & a_FCentroids) const
{
  a_rhsAtNodes.resize(m_grid.getNumInteriorNodes());
  
  //zero out all entries in RHS
  for (int i = 0; i < a_rhsAtNodes.size(); i ++)
    {
      a_rhsAtNodes[i] = 0.0;
    }

  //loop through each element
  for(int iEl = 0; iEl < m_grid.getNumElts(); iEl ++)
    {
      //access element with id iEL
      const Element& e = m_grid.element(iEl);

      //loop through the local vertices of element e
      for(int iVert = 0; iVert < VERTICES; iVert ++)
        {
          //xn is the corresponding vertex in element e
          const Node& xn = m_grid.node(e[iVert]);

          //we only populate b_n if xn is interior
          if (xn.isInterior())
            {
              //Using the pseudocode from page 22 of lecture 8. We know that Psi_n^h(centroid)=1/3 since it is piecewise linear
              a_rhsAtNodes[xn.getInteriorNodeID()] += m_grid.elementArea(iEl)*a_FCentroids[iEl]*(1/3);
            }
        }
    }

}


const FEGrid& FEPoissonOperator::getFEGrid() const
{
  return m_grid;
}


const SparseMatrix& FEPoissonOperator::matrix() const
{
  return m_matrix;
}
