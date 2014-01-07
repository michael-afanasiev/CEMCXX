#include <iostream>
#include "classes.h"
using namespace std;

int main () 
{
  
  int num_mesh_files;
  ifstream myfile;
  string * fname;
  
  Exodus_file exo;
  Mesh        msh;
  
  cout << "Begin model building.\n";  
  
  open_driver           ( myfile );
  read_driver_nMshFiles ( myfile, num_mesh_files );
  
  exo.openFile ();
  
  msh.getInfo      (exo.idexo);
  msh.allocateMesh (exo.idexo);
  msh.getCoord     (exo.idexo);
  
  exo.closeFile ();
  
  return 0;
  
}