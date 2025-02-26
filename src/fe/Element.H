#ifndef _ELEMENT_H_
#define _ELEMENT_H_
#define VERTICES 3
#include <array>
using namespace std;
/// A Class describing a triangular element in 2D. 
/** 
This class stores a an array of three integers that index into locations
in the vector containing the Nodes. The rest of the information defining the element (location in physical space, areas, centroids...) is computed from the node data.
 */
class Element
{
public: 
  Element();
  /// Constructor. 
  Element(array<int, VERTICES> a_vertices);

  /// local indexing to get nodeNumber. 
  /** Possible indices are 0,1,2, ordered clockwise around the center of the element.
   */
  const int& operator[](const int& a_localNodeNumber) const;
  /// Access to all the vertices at once.
  const array<int, VERTICES>& vertices() const;

private:
  array<int, VERTICES> m_vertices;
};
#endif
