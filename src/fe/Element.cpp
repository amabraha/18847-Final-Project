#include <cassert>
#include "Element.H"

Element::Element(array<int, VERTICES> a_vertices)
{
  m_vertices = a_vertices;
};


Element::Element()
{
 
};
const int& Element::operator[](const int& a_localNodeNumber) const
{
  assert(a_localNodeNumber < VERTICES);
  return m_vertices[a_localNodeNumber];
};


const array<int, VERTICES>& Element::vertices() const
{
  return m_vertices;
}
