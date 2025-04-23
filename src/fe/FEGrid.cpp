#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>
#include <vector>
#include "Node.H"
#include "Element.H"
#include "FEGrid.H"
#include "VisitWriter.H"
#include <fstream>

FEGrid::FEGrid() : m_numInteriorNodes(0)
{
}

FEGrid::FEGrid(const std::string &a_nodeFileName, const std::string &a_elementFileName)
{
  ifstream nodes(a_nodeFileName.c_str());
  int ncount, dim, attributes, boundaryMarkers;
  nodes >> ncount >> dim >> attributes >> boundaryMarkers;
  // cout<<ncount<<" "<<dim<<" "<<endl;

  m_nodes.resize(ncount);
  m_numInteriorNodes = 0;
  for (int i = 0; i < ncount; i++)
  {
    int vertex, type;
    array<double, DIM> x;
    nodes >> vertex >> x[0] >> x[1] >> type;
    vertex--;
    if (type == 1)
    {
      m_nodes[vertex] = Node(x, -1, false);
    }
    else
    {
      m_nodes[vertex] = Node(x, m_numInteriorNodes, true);
      m_numInteriorNodes++;
    }
  }

  ifstream elements(a_elementFileName.c_str());
  int ncell, nt;
  elements >> ncell >> nt >> attributes;
  array<int, VERTICES> vert;
  m_elements.resize(ncell);
  for (int i = 0; i < ncell; i++)
  {
    int cellID;
    elements >> cellID >> vert[0] >> vert[1] >> vert[2];
    vert[0]--;
    vert[1]--;
    vert[2]--;
    cellID--;
    m_elements[cellID] = Element(vert);
  }
};
array<double, DIM> FEGrid::gradient(
    const int &a_eltNumber,
    const int &a_nodeNumber) const
{
  const Element &e = m_elements[a_eltNumber];
  const Node &n = m_nodes[e[a_nodeNumber]];
  // assert(n.isInterior()); 
  struct xb
  {
    double x[DIM];
  };
  array<double, 2> xbase = n.getPosition();
  array<array<double, DIM>, VERTICES - 1> dx;
  for (int ivert = 0; ivert < VERTICES - 1; ivert++)
  {
    int otherNodeNumber = e[(a_nodeNumber + ivert + 1) % VERTICES];
    dx[ivert] = m_nodes[otherNodeNumber].getPosition();
    for (int idir = 0; idir < DIM; idir++)
    {
      dx[ivert][idir] -= xbase[idir];
    }
  }
  // WARNING: the following calculation is correct for triangles in 2D *only*.
  double det = dx[0][0] * dx[1][1] - dx[1][0] * dx[0][1];
  array<double, DIM> retval;

  retval[0] = (-(dx[1][1] - dx[0][1]) / det);
  retval[1] = (-(dx[1][0] - dx[0][0]) / det);
  return retval;
};
array<double, DIM> FEGrid::centroid(const int &a_eltNumber) const
{
  const Element &e = m_elements[a_eltNumber];
  array<double, DIM> retval;
  for (int i = 0; i < DIM; i++)
  {
    retval[i] = 0.0;
  }

  for (int ivert = 0; ivert < VERTICES; ivert++)
  {
    const Node &n = m_nodes[e[ivert]];
    array<double, DIM> x = n.getPosition();
    for (int idir = 0; idir < DIM; idir++)
    {
      retval[idir] += x[idir];
    }
  }
  for (int idir = 0; idir < DIM; idir++)
  {
    retval[idir] /= VERTICES;
  }
  return retval;
}
double FEGrid::elementArea(const int &a_eltNumber) const
{
  const Element &e = m_elements[a_eltNumber];
  const Node &n = m_nodes[e[0]];
  array<double, DIM> xbase = n.getPosition();
  array<array<double, VERTICES - 1>, DIM> dx;
  for (int ivert = 1; ivert < VERTICES; ivert++)
  {
    int otherNodeNumber = e[ivert];
    dx[ivert - 1] = m_nodes[otherNodeNumber].getPosition();
    for (int idir = 0; idir < DIM; idir++)
    {
      dx[ivert - 1][idir] -= xbase[idir];
    }
  }
  // WARNING: the following calculation is correct for triangles in 2D *only*.
  double area = fabs(dx[0][0] * dx[1][1] - dx[1][0] * dx[0][1]) / 2;
  return area;
}

// double FEGrid::elementValue(const int& a_eltNumber,
// 			   const int& a_localNodeNumber) const
// {
//   const Element& e = m_elements[a_eltNumber];
//   const Node& n=m_nodes[e[a_localNodeNumber]];
//   assert(n.isInterior());

//   double xbase[DIM];
//   n.getPosition(xbase);
//   double value = 0;

//   for (int idir = 0; idir < DIM; idir++)
//     {
//       value += a_xVal[idir] - xbase[idir];
//     }
//   return value;
// }
const Node &FEGrid::getNode(const int &a_eltNumber, const int &a_localNodeNumber) const
{
  return m_nodes[m_elements[a_eltNumber][a_localNodeNumber]];
}
int FEGrid::getNumElts() const
{
  return m_elements.size();
}

int FEGrid::getNumNodes() const
{
  return m_nodes.size();
}

int FEGrid::getNumInteriorNodes() const
{
  return m_numInteriorNodes;
}

const Element &FEGrid::element(int i) const
{
  return m_elements[i];
}
const Node &FEGrid::node(int i) const
{
  return m_nodes[i];
}

const char *FEWrite(FEGrid *a_grid, const char *a_filename)
{
  vector<double> data(a_grid->getNumNodes(), 0.0);
  return FEWrite(a_grid, &data, a_filename);
}

int fileCount = 0;
const char *FEWrite(FEGrid *a_grid)
{
  char filename[10];
  snprintf(filename, 10, "grid%d.vtk", fileCount);
  return FEWrite(a_grid, filename);
}

const char *FEWrite(FEGrid *a_grid, vector<double> *a_scalarField)
{
  char filename[10];
  snprintf(filename, 10, "FEData%d.vtk", fileCount);
  return FEWrite(a_grid, a_scalarField, filename);
}

const char *FEWrite(FEGrid *a_grid, vector<double> *a_scalarField, const char *a_filename)
{
  int nNodes = a_scalarField->size();
  assert(a_grid->getNumNodes() == nNodes);
  vector<float> scalarField;
  for (int k = 0; k < nNodes; k++)
  {
    scalarField.push_back((*a_scalarField)[k]);
  }
  float *vars[1];
  vars[0] = &(scalarField[0]);
  int vardim[1] = {1};
  int centering[1] = {1};
  const char *const varnames[] = {"nodeData"};

  vector<float> pts(3 * nNodes);
  for (int i = 0; i < nNodes; i++)
  {
    int p = 3 * i;
    array<double, DIM> x = a_grid->node(i).getPosition();
    pts[p] = x[0];
    pts[p + 1] = x[1];
    pts[p + 2] = 0.0;
  }

  int ncell = a_grid->getNumElts();
  vector<int> cellType(ncell, VISIT_TRIANGLE);
  vector<int> conns(3 * ncell);
  for (int i = 0; i < ncell; i++)
  {
    array<int, VERTICES> vertices = a_grid->element(i).vertices();
    int e = 3 * i;
    conns[e] = vertices[0];
    conns[e + 1] = vertices[1];
    conns[e + 2] = vertices[2];
  }

  write_unstructured_mesh(a_filename, 0, nNodes, &(pts[0]), ncell,
                          &(cellType[0]), &(conns[0]), 1, vardim, centering,
                          varnames, vars);
  return a_filename;
}
