#include <iostream>
#include <fstream>
#include <vector>
#include <exodusII.h>

class Constants
{
  
public:
  double PI      = 3.1415926535897932;
  double R_EARTH = 6371.0;
  double o80     = 180.0;
};


class Model_file
{
  
private:
  
  void readSES3D     ();
  void readSPECFEM3D ();
  
  void populateSES3D ( std::string name, int &num_regions, int &num_params, 
    std::vector<std::vector<double>> &vec , char ftype );
  
public:
  
  int * region;
  int num_regions;
  int num_x;
  int num_y;
  int num_z;
  int num_p;
  
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<std::vector<double>> col_rad;
  std::vector<std::vector<double>> lon_rad;
  std::vector<std::vector<double>> col_deg;
  std::vector<std::vector<double>> lon_deg;
  std::vector<std::vector<double>> rad;
  std::vector<std::vector<double>> vsh;
  std::vector<std::vector<double>> vsv;
  std::vector<std::vector<double>> rho;
  std::vector<std::vector<double>> vpp;
  
  std::string input_model_directory;
  std::string input_model_file_type;
  std::string input_model_physics;
  std::string absolute_or_perturb;
  
  // Internal functions.
  
  void read                 ();
  void colLonRad2xyzSES3D   ();
  void populateRadiansSES3D ();
  
};

class Driver
{
  
private:
  
  void getToken ( std::string &test );
  
public:
 
  int N_HEADER = 1;  

  
  std::string * params;
  
  void closeDriver ( std::ifstream &myfile );
  void openDriver  ( std::ifstream &myfile );
  void readDriver  ( std::ifstream &myfile );
  
};

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
  int num_mesh_files;
  
  float vers;
  
  std::string fname;
  
  // Internal functions.
  
  void openFile ( std::string fname );
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
  void allocateMesh ();
  void deallocateMesh ();
           
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