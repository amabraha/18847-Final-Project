#include <cstdio>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <vector>
#include <array>
#include "Node.H"
using namespace std;
Node::Node()
{
  for (int idir = 0;idir < DIM;idir++)
    {
      m_position[idir] = DBL_MAX;
    }
  m_isInterior = true;
  m_interiorNodeID = -1;
};
  
Node::Node(array<double, DIM> a_position,
           const int& a_interiorNodeID,
           const bool& a_isInterior)
{
      m_position = a_position;
      m_isInterior = a_isInterior;
      m_interiorNodeID = a_interiorNodeID;
};
const array<double, DIM>& Node::getPosition() const
{
  return m_position;
};

const int& Node::getInteriorNodeID() const
{
  return m_interiorNodeID;
};

const bool& Node::isInterior() const
{
  return m_isInterior;
};


  
