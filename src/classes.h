#include <exodusII.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

void open_driver           ( std::ifstream &myfile );
void read_driver_nMshFiles ( std::ifstream &myfile, int &num_mesh_files );  

class Exodus_file
{
public:  
  
  int relert;
  int chrret; 
  int ier; 
  int idexo; 
  int nump;
  int comp_ws;
  int io_ws;
  
  float vers;
  
  // Internal functions.
  
  void openFile ();
  void closeFile ();
  
};

class Mesh
{
public:
  
  int num_nodes;
  int num_elem;
  int ier;
  int num_dim;
  int num_elem_blk;
  int num_node_sets;
  int num_side_sets;
  
  double * c11;
  double * c12;
  double * c13;
  double * c14;
  double * c15;
  double * c16;
  double * c22;
  double * c23;
  double * c24;
  double * c25;
  double * c26;
  double * c33;
  double * c34;
  double * c35;
  double * c36;
  double * c44;
  double * c45;
  double * c46;
  double * c55;
  double * c56;
  double * c66;
  double * xmsh;
  double * ymsh;
  double * zmsh;
  
  char name;
  char title [MAX_LINE_LENGTH+1]; 
  
  // Internal functions.
  
  void getInfo (int exoid);
  void getCoord (int exoid);
  void allocateMesh (int exoid);
           
};

class Block 
{
  
  int n_elem_glob;
  
  double col_min;
  double col_max;
  double lon_min;
  double lon_max;
  double rad_min;
  double rad_max;
  
  double * var;
  
};