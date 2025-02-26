
#include "VisitWriter.H"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2)
    {
      cout << "this program takes one argument that is the .node and .ele ";
      cout << "file prefix"<<endl;
      return 1;
    }
  string prefix(argv[1]);
  string file = prefix+".node";
  ifstream nodes(file.c_str());
  int ncount, dim, attributes, boundaryMarkers;
  nodes>>ncount>>dim>>attributes>>boundaryMarkers;
  cout<<ncount<<" "<<dim<<" "<<endl;

  vector<float> pts(ncount*3,0.0);
  for(int i=0; i<ncount; i++)
    {
      int vertex, type;
      double x; double y;
      nodes>>vertex>>x>>y>>type;
      vertex--;
      pts[vertex*3] = x;
      pts[vertex*3+1] = y;
    }
  
  file = prefix+".ele";
  ifstream elements(file.c_str());
  int ncell, nt;
  elements>>ncell>>nt>>attributes;
  vector<int> conn(ncell*3);
  for(int i=0; i<ncell; i++)
    {
      int cellID, p1, p2, p3;
      elements>>cellID>>p1>>p2>>p3;
      cellID--;
      conn[cellID*3]=p1-1;
      conn[cellID*3+1] = p2-1;
      conn[cellID*3+2] = p3-1;
    }
  float* vars[1];
  vector<float> data(ncell,0.0);
  vars[0] = &(data[0]);
  int vardim[1] = {1};
  int centering[1] = {1};
  const char * const varnames[] = { "nodeData" };

  vector<int> cellType(ncell, VISIT_TRIANGLE);

  string fname = "ShowGridOutput.vtk";
  write_unstructured_mesh(fname.c_str(), 0, ncount, &(pts[0]), ncell, 
			  &(cellType[0]), &(conn[0]), 1, vardim, centering,
			  varnames, vars);
  return 0;
}
