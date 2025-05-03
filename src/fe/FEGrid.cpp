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
#include <iostream>

/* Triangle related */
#define REAL double
#define VOID void
extern "C" {
  #define ANSI_DECLARATORS
  #include "triangle.h"
}

FEGrid::FEGrid() : m_numInteriorNodes(0)
{
}

/* Helper for 3D extrusion */
int FEGrid::add_upper_node(const Node &node, double height)
{
  array<double, DIM> x;
  array<double, DIM> x_0 = node.getPosition();
  
  // Grow the layer
  x[0] = x_0[0];
  x[1] = x_0[1];
  x[2] = x_0[2] + height;

  int in_id = -1;
  bool is_in = node.isInterior();
  if (is_in) {
    in_id = m_numInteriorNodes;
    m_numInteriorNodes++;
  }

  Node node_new = Node(x, in_id, is_in);
  // Add new node
  m_nodes.push_back(node_new);
  return m_nodes.size() - 1;
}

FEGrid::FEGrid(const std::string &a_nodeFileName, const std::string &a_elementFileName)
{
  ifstream nodes(a_nodeFileName.c_str());
  int ncount, dim, attributes, boundaryMarkers;
  nodes >> ncount >> dim >> attributes >> boundaryMarkers;
  // cout<<ncount<<" "<<dim<<" "<<endl;

  // read nodes
  m_nodes.resize(ncount);
  m_numInteriorNodes = 0;
  for (int i = 0; i < ncount; i++)
  {
    int vertex, type;
    array<double, DIM> x;
    nodes >> vertex >> x[0] >> x[1] >> type;
    if (DIM == 3)
    {
      x[2] = 0.0;
    }
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
  // read elements
  ifstream elements(a_elementFileName.c_str());
  int ncell, nt;
  elements >> ncell >> nt >> attributes;
  array<int, VERTICES> vert;
  m_elements.resize(ncell);
  if (DIM == 3)
  {
    // Add new elements for 3d extrusion
    m_elements.resize(ncell * 3);
  }
  for (int i = 0; i < ncell; i++)
  {
    int cellID;
    elements >> cellID >> vert[0] >> vert[1] >> vert[2];
    vert[0]--;
    vert[1]--;
    vert[2]--;
    cellID--;

    // 2d extrude into 3d, grow a prism layer and split into 3 tetras
    if (DIM == 3) 
    {
      Node& v1 = m_nodes[vert[0]];
      Node& v2 = m_nodes[vert[1]];
      Node& v3 = m_nodes[vert[2]];
      array<double, DIM> x_1 = v1.getPosition();
      array<double, DIM> x_2 = v2.getPosition();
      array<double, DIM> x_3 = v3.getPosition();

      // Compute 2D area (cross product of vectors)
      double area = 0.5 * fabs((x_2[0]-x_1[0])*(x_3[1]-x_1[1])-(x_3[0]-x_1[0])*(x_2[1] - x_1[1]));

      // Compute height to preserve h: the 3D tetrahedron volumes have a cube root equivalent to the square root of the original 2D elements
      double h = 3.0 * sqrt(area);

      // Grow the upper layer into a prism
      int v4_id = add_upper_node(v1, h);
      int v5_id = add_upper_node(v2, h);
      int v6_id = add_upper_node(v3, h);
      // cout << "v4_id " << v4_id << endl;

      // Divide into 3 tetras
      // ref: https://www.researchgate.net/publication/221561839_How_to_Subdivide_Pyramids_Prisms_and_Hexahedra_into_Tetrahedra
      // Tetra 1: v1, v4, v5, v6
      // Tetra 2: v1, v2, v5, v6
      // Tetra 3: v1, v2, v3, v6
      array<int, VERTICES> tet1;
      tet1[0] = vert[0];
      tet1[1] = v4_id;
      tet1[2] = v5_id;
      tet1[3] = v6_id;
      array<int, VERTICES> tet2;
      tet2[0] = vert[0];
      tet2[1] = vert[1];
      tet2[2] = v5_id;
      tet2[3] = v6_id;
      array<int, VERTICES> tet3;
      tet3[0] = vert[0];
      tet3[1] = vert[1];
      tet3[2] = vert[2];
      tet3[3] = v6_id;
      

      // Add all tetras
      m_elements[cellID * 3] = Element(tet1);
      m_elements[cellID * 3 + 1] = Element(tet2);
      m_elements[cellID * 3 + 2] = Element(tet3);
    }
    else {
      m_elements[cellID] = Element(vert);
    }
  }
};

/**
  Constructor to Integrate Triangle library
  1. Read the .poly file for nodes, segments
  2. call triangle to generate nodes and elements
  3. initialize nodes and elements
*/
FEGrid::FEGrid(const std::string& a_polyFileName, const double max_area)
{
  struct triangulateio in, out;
  in.numberofcorners = 3;

  ifstream poly(a_polyFileName.c_str());

  if (!poly) 
  {
      std::cout << "Error: File does not exist or failed to open: " << a_polyFileName << std::endl;
      exit(1);
  }

  // Read the vertices: build up pointlist
  // header of poly file
  int ncount, dim, attributes, boundaryMarkers;
  poly >> ncount >> dim >> attributes >> boundaryMarkers;
  // initialize data structures
  in.numberofpoints = ncount;
  in.numberofpointattributes = attributes;
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
  if (in.numberofpointattributes > 0) {
    in.pointattributelist = (REAL *) malloc(in.numberofpoints *
                                          in.numberofpointattributes *
                                          sizeof(REAL));
  } else
  {
    in.pointattributelist = (REAL *) NULL;
  }
  if (boundaryMarkers > 0) {
    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
  } else
  {
    in.pointmarkerlist = (int *) NULL;
  }

  // body
  for(int i=0; i<ncount; i++)
    {
      int vertex;
      double x, y;
      poly >> vertex >> x >> y;
      vertex--;
      // read position
      in.pointlist[vertex * 2] = x;
      in.pointlist[vertex * 2 + 1] = y;

      // read attributes
      for (int j = 0; j < in.numberofpointattributes; j++)
      {
        poly >> in.pointattributelist[vertex*in.numberofpointattributes + j];
      }

      // read boundary marker
      if (boundaryMarkers > 0) {
        poly >> in.pointmarkerlist[vertex];
      }
    }
  

  in.trianglelist = (int *) NULL;
  
  // Read the segments: build up segment list
  int segcount, segboundaryMarkers;
  // header
  poly >> segcount >> segboundaryMarkers;
  // initialize data structures
  in.numberofsegments = segcount;
  in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  if (segboundaryMarkers > 0) {
    in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
  } else
  {
    in.segmentmarkerlist = (int *) NULL;
  }
  // body
  for(int i=0; i<segcount; i++)
    {
      int segID, endpoint1, endpoint2;
      poly >> segID >> endpoint1 >> endpoint2;
      segID--;
      // read boundary marker
      if (segboundaryMarkers > 0) {
        poly >> in.segmentmarkerlist[segID];
      }
      // read segment endpoints (node ID)
      in.segmentlist[segID * 2] = endpoint1;
      in.segmentlist[segID * 2 + 1] = endpoint2;
    }
  
  // Read the holes
  int holecount;
  // header
  poly >> holecount;
  in.numberofholes = holecount;
  // body
  if (in.numberofholes > 0)
  {
    in.holelist = (REAL *)malloc(in.numberofholes * 2 * sizeof(REAL));
    for(int i = 0; i < holecount; i++)
      {
        int holeID;
        REAL x, y;
        poly >> holeID >> x >> y;
        holeID--;
        in.holelist[holeID * 2] = x;
        in.holelist[holeID * 2 + 1] = y;
      }
  }
  else {
    in.holelist = (REAL *) NULL; 
  }

  // Assumed fixed max area, no regional attributes and/or area constraints in .poly
  in.numberofregions = 0;
  in.regionlist = (REAL *) NULL;

  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out'.                                    */

  out.pointlist = (REAL *) NULL; 
  out.pointattributelist = (REAL *) NULL;
  out.trianglelist = (int *) NULL;
  out.triangleattributelist = (REAL *) NULL;
  out.pointmarkerlist = (int *) NULL;
  out.segmentlist = (int *) NULL;
  out.segmentmarkerlist = (int *) NULL;

  /* Refine the triangulation according to the attached */
  /*   triangle area constraints.                       */
  char switches[20];
  snprintf(switches, sizeof(switches), "pqa%f", max_area);
  triangulate(switches, &in, &out, (struct triangulateio *) NULL);

  // Read nodes from pointlist
  ncount = out.numberofpoints;
  boundaryMarkers = (out.pointmarkerlist != NULL);
  m_nodes.resize(ncount);
  m_numInteriorNodes= 0;
  for(int vertex=0; vertex<ncount; vertex++)
    {
      int type = out.pointmarkerlist[vertex];
      array<double, DIM> x;
      x[0] = out.pointlist[vertex * 2];
      x[1] = out.pointlist[vertex * 2 + 1];
      if (DIM == 3)
      {
        x[2] = 0.0;
      }
      if(type == 1)
      {
        m_nodes[vertex] = Node(x,-1, false);
      }
      else
      {
        m_nodes[vertex] = Node(x, m_numInteriorNodes, true);
        m_numInteriorNodes++;
      }
    }

  // Read elements from trianglelist
  int ncell;
  ncell = out.numberoftriangles;
  m_elements.resize(ncell);
  if (DIM == 3)
  {
    // Add new elements for 3d extrusion
    m_elements.resize(ncell * 3);
  }
  for(int i = 0; i < ncell; i++)
    {
      int cellID = i;
      array<int, VERTICES> vert;
      // read vertices (node ID)
      vert[0] = out.trianglelist[cellID * 3] - 1; // 1-based indexing
      vert[1] = out.trianglelist[cellID * 3 + 1] - 1;
      vert[2] = out.trianglelist[cellID * 3 + 2] - 1;
      // 2d extrude into 3d, grow a prism layer and split into 3 tetras
      if (DIM == 3) 
      {
        Node v1 = m_nodes[vert[0]];
        Node v2 = m_nodes[vert[1]];
        Node v3 = m_nodes[vert[2]];
        array<double, DIM> x_1 = v1.getPosition();
        array<double, DIM> x_2 = v2.getPosition();
        array<double, DIM> x_3 = v3.getPosition();

        // Compute 2D area (cross product of vectors)
        double area = 0.5 * fabs((x_2[0]-x_1[0])*(x_3[1]-x_1[1])-(x_3[0]-x_1[0])*(x_2[1] - x_1[1]));

        // Compute height to preserve h: the 3D tetrahedron volumes have a cube root equivalent to the square root of the original 2D elements
        double h = 3.0 * sqrt(area);

        // Grow the upper layer into a prism
        int v4_id = add_upper_node(v1, h);
        int v5_id = add_upper_node(v2, h);
        int v6_id = add_upper_node(v3, h);
        // cout << "v6 id " << v6_id << endl;

        // Divide into 3 tetras
        // ref: https://www.researchgate.net/publication/221561839_How_to_Subdivide_Pyramids_Prisms_and_Hexahedra_into_Tetrahedra
        // Tetra 1: v1, v4, v5, v6
        // Tetra 2: v1, v2, v5, v6
        // Tetra 3: v1, v2, v3, v6
        array<int, VERTICES> tet1;
        tet1[0] = vert[0];
        tet1[1] = v4_id;
        tet1[2] = v5_id;
        tet1[3] = v6_id;
        array<int, VERTICES> tet2;
        tet2[0] = vert[0];
        tet2[1] = vert[1];
        tet2[2] = v5_id;
        tet2[3] = v6_id;
        array<int, VERTICES> tet3;
        tet3[0] = vert[0];
        tet3[1] = vert[1];
        tet3[2] = vert[2];
        tet3[3] = v6_id;

        // Add all tetras
        m_elements[cellID * 3] = Element(tet1);
        m_elements[cellID * 3 + 1] = Element(tet2);
        m_elements[cellID * 3 + 2] = Element(tet3);
      }
      else {
        m_elements[cellID] = Element(vert);
      }
    }

  // Clean up
  free(in.pointlist);
  free(in.pointattributelist);
  if (boundaryMarkers > 0) {
    free(in.pointmarkerlist);
  }
  if (in.numberofsegments > 0) {
    free(in.segmentlist);
    if (segboundaryMarkers > 0) {
      free(in.segmentmarkerlist);
    }
  }
  if (in.numberofholes > 0) {
    free(in.holelist);
  }
  
  free(out.pointlist);
  free(out.pointattributelist);
  free(out.trianglelist);
  free(out.triangleattributelist);
};

/** Gradient of linear basis function
*/
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
  array<double, DIM> xbase = n.getPosition();
  array< array<double, DIM> , VERTICES - 1> dx;
  for (int ivert = 0;ivert < VERTICES - 1; ivert++)
    {
      // Build vectors to other nodes of the element
      int otherNodeNumber = e[(a_nodeNumber + ivert + 1)%VERTICES];
      dx[ivert] = m_nodes[otherNodeNumber].getPosition();
      for (int idir = 0;idir < DIM;idir++)
        {
          dx[ivert][idir] -=xbase[idir];
        }
    }        
  if (DIM == 2) {
    // WARNING: the following calculation is correct for triangles in 2D *only*.
    // Determinant (for calculating matrix inverse, 2*area)
    double det = dx[0][0] * dx[1][1] - dx[1][0] * dx[0][1];
    array<double, DIM> retval;
    // Solve for gradient
    retval[0] = (-(dx[1][1] - dx[0][1]) / det);
    retval[1] = (-(dx[1][0] - dx[0][0]) / det);
    return retval;
  } 
  else {
    // Tetrahedrons in 3D
    assert(DIM==3);
    // Determinant (dot product, cross product, 6*volume)
    double det = dx[0][0] * (dx[1][1] * dx[2][2] - dx[1][2] * dx[2][1]) 
               - dx[0][1] * (dx[1][0] * dx[2][2] - dx[1][2] * dx[2][0])
               + dx[0][2] * (dx[1][0] * dx[2][1] - dx[1][1] * dx[2][0]);
    // Solve for gradient in 3D (compute inverse, multiply by -1)
    array<double, DIM> retval;
    retval[0] = (-(dx[1][1] * dx[2][2] - dx[1][2] * dx[2][1]) / det);
    retval[1] = (-(dx[1][0] * dx[2][2] - dx[1][2] * dx[2][0]) / det);
    retval[2] = (-(dx[1][0] * dx[2][1] - dx[1][1] * dx[2][0]) / det);
    return retval;
  }
};

array<double, DIM> FEGrid::centroid(const int& a_eltNumber) const
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


double FEGrid::elementArea(const int& a_eltNumber) const
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
      int otherNodeNumber = e[ivert];
     dx[ivert-1] = m_nodes[otherNodeNumber].getPosition();
      for (int idir = 0;idir < DIM;idir++)
        {
          dx[ivert-1][idir] -=xbase[idir];
        }
    }  
  }      
  if (DIM == 2) 
  {
    // WARNING: the following calculation is correct for triangles in 2D *only*.
    double area = fabs(dx[0][0]*dx[1][1] - dx[1][0]*dx[0][1])/2.0;
    return area;
  } else 
  {
    assert(DIM == 3);
    // Calculate volume for 3D tetrahedron (absolute value of Determinant / 6)
    double volume = fabs(dx[0][0]*(dx[1][1]*dx[2][2]-dx[1][2]*dx[2][1]) 
                - dx[0][1]*(dx[1][0]*dx[2][2]-dx[1][2]*dx[2][0])
                + dx[0][2]*(dx[1][0]*dx[2][1]-dx[1][1]*dx[2][0]))/6.0;
    return volume;
  }
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
  const char * const varnames[] = { "nodeData" };
  
  // Always 3d coordinates
  vector<float> pts(3 * nNodes);
  for(int i = 0; i < nNodes; i++)
    {
      int p = 3 * i;
      array<double, DIM> x = a_grid->node(i).getPosition();
      pts[p] = x[0];
      pts[p + 1] = x[1];
      if (DIM == 2) 
      {
        pts[p + 2] = 0.0;
      }
      else {
        pts[p + 2] = x[2];
      }
    }

  int ncell = a_grid->getNumElts();
  vector<int> cellType;
  if (DIM == 2) {
    // 2D: 3 points per element
    cellType.resize(ncell, VISIT_TRIANGLE);
  } else {
    // 3D: 4 points per element
    cellType.resize(ncell, VISIT_TETRA);
  }

  vector<int> conns(VERTICES*ncell);
  for(int i = 0; i < ncell; i++)
    {
      array<int, VERTICES> vertices = a_grid->element(i).vertices();
      int e = VERTICES * i;
      for (int j = 0; j < VERTICES; j++) {
        conns[e + j] = vertices[j];
      }
    }

  write_unstructured_mesh(a_filename, 0, nNodes, &(pts[0]), ncell, 
			  &(cellType[0]), &(conns[0]), 1, vardim, centering,
			  varnames, vars);
  return a_filename;
}
