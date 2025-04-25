#ifndef _FEGRID_H_
#define _FEGRID_H_

#include <cstdio>   
#include <cmath>    
#include <cassert>
#include <cstring> 
#include <vector> 
#include <string>
#include "Node.H"   
#include "Element.H"

using namespace std;

class FEGrid
{
public:
  FEGrid();
  /// Constructor by reading from file.
  FEGrid(const std::string& nodeFile, const std::string& a_elementFileName);
  FEGrid(const std::string& a_polyFileName, const double max_area);

  /// Compute gradient of basis function at node a_nodeNumber = 0,..,DIM on element a_eltNumber.
  array<double, DIM> gradient(
                             const int& a_eltNumber,
                             const int& a_localNodeNumber) const;
  /// Compute centroid of element.                                
  array<double, DIM> centroid(const int& a_eltNumber) const;
  /// Compute area of element.
  double elementArea(const int& a_eltNumber) const;
 
  /// get reference to node, given an element number and a local node number.                      
  const Node& getNode(const int& a_eltNumber,
                      const int& a_localNodeNumber) const; 
  /// Get number of elements.
  int getNumElts() const;
  int getNumNodes() const;
  int getNumInteriorNodes() const;
  /// Access an element by indexing into  m_elements.
  const Element& element(int i) const;
  /// Access a node by indexing into m_nodes.
  const Node& node(int i) const;
  int add_upper_node(const Node &node, double height);
  private:                            
  vector<Node > m_nodes;
  vector<Element > m_elements;
  int m_numInteriorNodes;
};

const char* FEWrite(FEGrid* a_grid, const char* a_filename);
const char* FEWrite(FEGrid* a_grid);
const char* FEWrite(FEGrid* a_grid, vector<double>* a_scalarField, const char* a_filename);
const char* FEWrite(FEGrid* a_grid, vector<double>* a_scalarField);


#endif
