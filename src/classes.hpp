#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <exodusII.h>

#include "kdtree.h"

class Mesh;
class Driver;
class Exodus_file;
class Interpolator;

class Region
{

public:
  std::vector < Exodus_file > regionsExo;

};

class Constants
{
  
public:
  
  double PI            = 3.141592653589793;
  double PIo2          = 1.570796326794896;
  double R_EARTH       = 6371.0;
  double o80           = 180.0;
  double ninty         = 90.0;
  double aniCorrection = 0.188078;
  
};

class Model_file
{
  
public:
  
  void readSES3D     ();
  void readSPECFEM3D ();
  
  void populateSES3D ( std::string name, int &num_regions, int &num_params, 
    std::vector<std::vector<double>> &vec , char ftype );
      
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
  
  double colMin = 180.;
  double colMax = 0.;
  double lonMin = 360.;
  double lonMax = -360.;
  double radMin = 6371.;
  double radMax = 0.;
  
  std::vector <double> x;
  std::vector <double> y;
  std::vector <double> z;
  
  double xMin;
  double xMax;
  double yMin;
  double yMax;
  double zMin;
  double zMax;

  double *c11;
  double *c12;
  double *c13;
  double *c14;
  double *c15;
  double *c16;
  double *c22;
  double *c23;
  double *c24;
  double *c25;
  double *c26;
  double *c33;
  double *c34;
  double *c35;
  double *c36;
  double *c44;
  double *c45;
  double *c46;
  double *c55;
  double *c56;
  double *c66;
  double *rhoMsh;
  
  double forRot11;
  double forRot21;
  double forRot31;
  double forRot12;
  double forRot22;
  double forRot32;
  double forRot13;
  double forRot23;
  double forRot33;

  double backRot11;
  double backRot21;
  double backRot31;
  double backRot12;
  double backRot22;
  double backRot32;
  double backRot13;
  double backRot23;
  double backRot33;
  
  bool colReg1    = false;
  bool colReg2    = false;
  bool lonReg1    = false;
  bool lonReg2    = false;
  bool lonReg3    = false;
  bool lonReg4    = false;
  bool wrapAround = false;
  bool doRotate   = false;
  
  kdtree *tree;

  std::vector <std::vector <double> > col_rad;
  std::vector <std::vector <double> > lon_rad;
  std::vector <std::vector <double> > col_deg;
  std::vector <std::vector <double> > lon_deg;
  std::vector <std::vector <double> > rad;
  std::vector <std::vector <double> > vsh;
  std::vector <std::vector <double> > vsv;
  std::vector <std::vector <double> > rho;
  std::vector <std::vector <double> > vpp;  
  
  std::string input_model_directory;
  std::string input_model_file_type;
  std::string input_model_physics;
  std::string absolute_or_perturb;
  std::string intentions;
  std::string dimensions;
  std::string interpolation;
  std::string output_model_physics;
  std::string mesh_directory;
  
  // Internal functions.
  
  void read                 ();
  void colLonRad2xyzSES3D   ();
  void populateRadiansSES3D ();
  void openUp               ();
  void populateParams       ( Driver &drv, Exodus_file &exo );
  void dePopulateSES3D      ( std::string, std::vector<std::vector<double>> );
  void writeSES3D           ();
  void populateRadians      ( std::vector < std::vector <double> > &deg, 
                              std::vector < std::vector <double> > &rad );
  void createKDTreeUnpacked ();
  void projectSubspace      ();
  
};

class Discontinuity
{
  
public:
  
  bool inCrust;
  
  std::vector <std::vector <double> > crust_col_deg;
  std::vector <std::vector <double> > crust_lon_deg;
  std::vector <std::vector <double> > crust_col_rad;
  std::vector <std::vector <double> > crust_lon_rad;
  std::vector <std::vector <double> > crust_vs;
  std::vector <std::vector <double> > crust_dp;
  
  kdtree *tree;
  
  void read               ();
  void createKDTreePacked ();
  void lookCrust          ( Mesh &msh, double &mshCol, double &mshLon, 
                            double &mshRad, int &mshInd );

};

class Utilities
{
  
public:
    
  double col2Lat        ( double &in, char flag );
  void rotateForward    ( double &x, double &y, double &z, double &xRot, 
                          double &yRot, double &zRot, Model_file &mod );   
  void rotateBackward   ( double &x, double &y, double &z, double &xRot, 
                          double &yRot, double &zRot, Model_file &mod );   
  void xyz2ColLonRadDeg ( double &x,   double &y,   double &z, 
                          double &col, double &lon, double &rad );
  void xyz2ColLonRadRad ( double &x,   double &y,   double &z, 
                          double &col, double &lon, double &rad );
  void colLonRadDeg2xyz ( double  col, double  lon, double  rad,
                          double &x,   double &y,   double &z );  
  void inquireRotate    ( Model_file &mod );
  void colLonRadRad2xyz ( double col,  double lon,  double rad,
                          double &x,   double &y,   double &z );                          
  void convertBary      ( double &xp, double &yp, double &zp, 
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
  int num_mesh_files;
  
  float vers;
  
  std::string fname;
  
  std::vector < std::string > colReg;
  std::vector < std::string > lonReg;
  std::vector < std::string > radReg;    
  
  std::vector < int > totalBlocks;
  
  // Internal functions.
  
  void openFile    ( std::string fname );
  void closeFile   ();
  void writeParams ( Mesh &msh );
  void merge       ( Region &reg, Model_file &mod );
  void splitBack   ();
  void writeSize   ( Mesh &msh );
  
};

class Mesh
{
public:
    
  int num_nodes;
  int num_elem;
  int ier;
  int num_dim;
  int num_elem_blk;
  int num_elem_in_blk;
  int num_node_sets;
  int num_side_sets;
  int num_node_per_elem=4;
  int num_nodes_in_elem;
  int num_attr;
  int exoid;
  
  double *c11;
  double *c12;
  double *c13;
  double *c14;
  double *c15;
  double *c16;
  double *c22;
  double *c23;
  double *c24;
  double *c25;
  double *c26;
  double *c33;
  double *c34;
  double *c35;
  double *c36;
  double *c44;
  double *c45;
  double *c46;
  double *c55;
  double *c56;
  double *c66;
  double *rho;
  double *Q__;
  double *elv;
  double *du1;
  double *du2;
  double *du3;
  
  double *xmsh;
  double *ymsh;
  double *zmsh;  
  double *siz;
  
  double colMin = 180.;
  double colMax = 0.;
  double lonMin = 180.;
  double lonMax = -180.;
  double radMin = 6371.;
  double radMax = 0.;
  
  kdtree *tree;
  
  std::multimap <int, std::vector <int> > elemOrder;  
  
  char name;
  char title [MAX_LINE_LENGTH+1]; 
  char elem_type [MAX_LINE_LENGTH+1];
  
  // Internal functions.
    
  void getInfo              ( int exoid, char mode );
  void populateCoord        ( );
  void populateParams       ( );
  void allocateMesh         ( );
  void reNormalize          ( Model_file &mod );
  void getConnectivity      ( int exoid );
  void deallocateMesh       ( Model_file &mod );
  void createKDTreeUnpacked ( );
  void getMinMaxRad         ( );
           
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

class Interpolator
{
  
public: 
    
  void interpolate  ( Mesh &msh, Model_file &mod, Discontinuity &dis );  
  double taper      ( double &x, double &y, double &z, Model_file &mod );
  void exterpolator ( Mesh &msh, Exodus_file &exo, Model_file &mod );
  int  recover      ( double &testX, double &testY, double &testZ, kdtree *tree,
                      Mesh &msh,
                      double &c11, double &c12, double &c13, double &c14, 
                      double &c15, double &c16, double &c22, double &c23, 
                      double &c24, double &c25, double &c26, double &c33, 
                      double &c34, double &c35, double &c36, double &c44, 
                      double &c45, double &c46, double &c55, double &c56, 
                      double &c66, double &rho, char mode ); 

};