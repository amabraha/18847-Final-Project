#include <cmath>  
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>
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
  /// Build the sparse matrix L
  m_matrix = SparseMatrix(m_grid.getNumInteriorNodes(), m_grid.getNumInteriorNodes());
  /// Iterate interior node pairs in elements
  for (int ele_id = 0; ele_id < m_grid.getNumElts(); ele_id++)
  {
    for (int local_node_i = 0; local_node_i < VERTICES; local_node_i++)
    {
      Node node_i = m_grid.getNode(ele_id, local_node_i);
      /// Check is interior
      if (node_i.isInterior())
      {
        array<double, DIM> gradient_i = m_grid.gradient(ele_id, local_node_i);
        for (int local_node_j = 0; local_node_j < VERTICES; local_node_j++)
        {
          Node node_j = m_grid.getNode(ele_id, local_node_j);
          if (node_j.isInterior())
          {
            array<double, DIM> gradient_j = m_grid.gradient(ele_id, local_node_j);
            double inner_product = 0;
            for (int dim = 0; dim < DIM; dim++)
            {
              inner_product += gradient_i[dim] * gradient_j[dim];
            }
            /// Fill up matrix incrementally
            array<int, 2> index = {node_i.getInteriorNodeID(), node_j.getInteriorNodeID()};
            m_matrix[index] += inner_product * m_grid.elementArea(ele_id);
          }
        }
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
              a_rhsAtNodes[xn.getInteriorNodeID()] += m_grid.elementArea(iEl)*a_FCentroids[iEl]*(1.0/3.0);
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
