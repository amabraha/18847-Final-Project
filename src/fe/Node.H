#ifndef _NODE_H_
#define _NODE_H_
#define DIM 2
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>
#include <array>
using namespace std;
class Node
{
 public:
  /// default constructor. 
  Node();
  /// Real constructor: location in space and flag indicating whether node is interior
  /// or on boundary.
  Node(array<double, DIM> a_position,
       const int&  a_interiorNodeID, 
       const bool& a_isInterior);
  /// Constant access to node Location in space.
  const array<double, DIM>&  getPosition() const;
  /// Index into interior node vector, provided the node is an interior node. Otherwise returns -1.
  const int& getInteriorNodeID() const;
  /// True if the node is an interior node.
  const bool& isInterior() const;

 private:
  array<double, DIM> m_position;
  bool  m_isInterior;
  int m_interiorNodeID;
};
#endif
