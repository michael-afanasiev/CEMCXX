#include <iostream>
#include <fstream>
#include <vector>
#include <exodusII.h>

class Mesh;

class Constants
{
  
public:
  
  double PI      = 3.141592653589793;
  double PIo2    = 1.570796326794896;
  double R_EARTH = 6371.0;
  double o80     = 180.0;
  double ninty   = 90.0;
  
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
  
  double rotAng;
  double rotRad;
  double rotVecX;
  double rotVecY;
  double rotVecZ;
  
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
  std::string intentions;
  std::string dimensions;
  std::string interpolation;
  
  // Internal functions.
  
  void read                 ();
  void colLonRad2xyzSES3D   ();
  void populateRadiansSES3D ();
  void openUp               ( Mesh &msh );
  
};

class Utilities
{
  
public:
  
  double col2Lat   ( double &in, char flag );
  void rotate      ( Model_file &mod );   
  void convertBary ( double &xp, double &yp, double &zp, 
                   double &x1, double &x2, double &x3, double &x4,
                   double &y1, double &y2, double &y3, double &y4,
                   double &z1, double &z2, double &z3, double &z4,
                   double &l1, double &l2, double &l3, double &l4 ); 
  
};


class Driver
{
  
private:
  
  void getToken ( std::string &test );
  
public:
 
  int N_HEADER = 1;  
  
  std::string *params;
  
  void closeDriver ( std::ifstream &myfile );
  void openDriver  ( std::ifstream &myfile );
  void readDriver  ( std::ifstream &myfile );
  
};

class Exodus_file
{
public:  

  int comp_ws = 8;
  int io_ws   = 0;
  
  int relert;
  int chrret; 
  int ier; 
  int idexo; 
  int nump;
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
  
  void getInfo          ( int exoid );
  void populateCoord    ( int exoid );
  void populateParams   ( int eoxid );
  void allocateMesh     ( int &num_nodes );
  void interpolateModel ( Model_file &mod );
  void deallocateMesh   ();
           
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