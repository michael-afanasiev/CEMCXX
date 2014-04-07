#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <exodusII.h>

#include "kdtree.h"

class Mesh;
class Driver;
class Region;
class Exodus_file;
class Interpolator;
class Mod1d;

class Region
{

public:
  std::vector < Exodus_file > regionsExo;

};

class Mod1d
{
  
public:
  
  void eumod ( double &, double &, double &, double & );
  
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
  double toMB          = 9.5367e-7;
  double tiny          = 1e-4;
  
  double innerCoreRad = 1221.0;
  double outerCoreRad = 3480.0;
  double R670         = 5701.0;
  double R400         = 5971.0; 
  
};

class Model_file
{
  
public:
  
  void readSES3D     ();
  void readSPECFEM3D ();
  
  void populateSES3D ( std::string name, int &num_regions, 
    std::vector <std::vector <double> > &vec , char ftype );
      
  int *region;
  int *KDdat;
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
  
  double refineSize;
  
  double colMin = 180.;
  double colMax = 0.;
  double lonMin = 360.;
  double lonMax = -360.;
  double radMin = 6371.;
  double radMax = 0.;
  
  std::vector <double> x;
  std::vector <double> y;
  std::vector <double> z;

  std::vector <double> c11;
  std::vector <double> c12;
  std::vector <double> c13;
  std::vector <double> c14;
  std::vector <double> c15;
  std::vector <double> c16;
  std::vector <double> c22;
  std::vector <double> c23;
  std::vector <double> c24;
  std::vector <double> c25;
  std::vector <double> c26;
  std::vector <double> c33;
  std::vector <double> c34;
  std::vector <double> c35;
  std::vector <double> c36;
  std::vector <double> c44;
  std::vector <double> c45;
  std::vector <double> c46;
  std::vector <double> c55;
  std::vector <double> c56;
  std::vector <double> c66;
  std::vector <double> rhoUnwrap;
  std::vector <double> vshUnwrap;
  std::vector <double> vsvUnwrap;
  std::vector <double> vppUnwrap;
  std::vector <double> radUnwrap;
      
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
  bool overwriteCrust;
  bool kernel;
  
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
  std::string intentions;
  std::string output_model_physics;
  std::string mesh_directory;
  
  // Internal functions.
  
  void read                 ( );
  void colLonRad2xyzSES3D   ( );
  void populateRadiansSES3D ( );
  void openUp               ( );
  void createKDTreeUnpacked ( Mesh &msh );
  void projectSubspace      ( );
  void writeSES3D           ( );
  void deallocate           ( );
  int  writeNetCDF          ( std::vector <std::vector <double>> &par,
                              std::string name );
  void populateParams       ( Driver &drv, Exodus_file &exo );
  void dePopulateSES3D      ( std::string, std::vector<std::vector<double>> );
  void populateRadians      ( std::vector < std::vector <double> > &deg, 
                              std::vector < std::vector <double> > &rad );
  
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
  
  std::vector <double> colElv;
  std::vector <double> lonElv;
  std::vector <double> elv;
  
  kdtree *crustTree;
  kdtree *elvTree;
  
  int *KDdatCrust;
  int *KDdatElv;
  
  void read                 ( );
  void createKDTreePacked   ( );
  void createKDTreeUnpacked ( );
  void deallocate           ( );
  void readTopography       ( );
  void lookCrust            ( Mesh &msh, double &mshCol, double &mshLon, 
                            double &mshRad, int &mshInd, bool &checkCrust );
  void lookTopo             ( Mesh &msh, double &mshCol, double &mshLon, 
                              double &mshRad, int &mshInd );

};

class Utilities
{
  
public:

  void inquireRotate    ( Model_file &mod );    
  double col2Lat        ( double &in, char flag );
  void colLonRadDeg2xyz ( double  col, double  lon, double  rad,
                          double &x,   double &y,   double &z );  
  void colLonRadRad2xyz ( double col,  double lon,  double rad,
                          double &x,   double &y,   double &z );                          
  void xyz2ColLonRadDeg ( double &x,   double &y,   double &z, 
                          double &col, double &lon, double &rad );
  void xyz2ColLonRadRad ( double &x,   double &y,   double &z, 
                          double &col, double &lon, double &rad );
  void rotateForward    ( double &x, double &y, double &z, double &xRot, 
                          double &yRot, double &zRot, Model_file &mod );   
  void rotateBackward   ( double &x, double &y, double &z, double &xRot, 
                          double &yRot, double &zRot, Model_file &mod );   
  void convertBary      ( double &xp, double &yp, double &zp, 
                          double &x1, double &x2, double &x3, double &x4,
                          double &y1, double &y2, double &y3, double &y4,
                          double &z1, double &z2, double &z3, double &z4,                          
                          double &l1, double &l2, double &l3, double &l4 );   
  void checkRegion      ( Mesh &msh, double &rad );
                          
  int getFilesize       ( std::string fname );
};


class Driver
{
  
  friend class Model_file;
    
private:
  
  int N_HEADER = 1;    
  std::string *params;
  
  void getToken    ( std::string &test );
  void closeDriver ( std::ifstream &myfile );
  
public: 
  
  void readDriver  ( std::ifstream &myfile );
  void checkUsage  ( Model_file &mod, std::string mode );
  void initialize  ( Model_file &mod, Discontinuity &dis, Utilities &utl, 
  Exodus_file &exo, Region &reg );
  void populateParams ( Model_file &mod );
  void report      ( Model_file &mod );
  
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
  
  bool allFiles = false;
  
  std::string fname;
  
  std::vector < std::string > colReg;
  std::vector < std::string > lonReg;
  std::vector < std::string > radReg;    
  
  std::vector < int > totalBlocks;
  
  // Internal functions.
  
  void writeNew   ( Mesh &msh );
  void merge       ( Region &reg, Model_file &mod );
  void openFile    ( std::string fname );
  void writeParams ( Mesh &msh );
  void closeFile   ( );
  
};

class Mesh
{
public:
    
  int ier;
  int exoid;
  int num_dim;
  int num_elem;
  int num_attr;
  int num_nodes;
  int num_elem_blk;
  int num_node_sets;
  int num_side_sets;
  int num_elem_in_blk;
  int num_nodes_in_elem;
  int numFound = 0;
  int num_node_per_elem=4;

  int *node_num_map;
  int *elem_num_map;
  int *KDdat;
  
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
  double *siz;
  
  double *xmsh;
  double *ymsh;
  double *zmsh;  
  
  double colMin = 180.;
  double colMax = 0.;
  double lonMin = 180.;
  double lonMax = -180.;
  double radMin = 6371.;
  double radMax = 0.;
  
  kdtree *tree;
  
  std::multimap <int, std::vector <int> > elemOrder;  
  std::vector < std::vector <int> > refineElemConn;
  
  char name;
  char title [MAX_LINE_LENGTH+1]; 
  char elem_type [MAX_LINE_LENGTH+1];
  
  std::string regString;
  
  // Internal functions.
    
  void getInfo                ( int exoid, char mode );
  void populateCoord          ( );
  void populateParams         ( );
  void allocateMesh           ( );
  void getConnectivity        ( int exoid );
  void deallocateMesh         ( Model_file &mod );
  void createKDTreeUnpacked   ( );
  void getMinMaxRad           ( );
  void getNodeNumMap          ( int exoid );
  void getElemNumMap          ( int exoid );
  void getElementConnectivity ( int exoid, Model_file &mod );  
           
};

class Interpolator
{
  
public: 
  
  std::vector <int> refineArr;
  std::vector < std::vector <int> > elemWithin;
    
  void findNodes        ( Mesh &msh, Model_file &mod, std::ofstream &myfile );
  void interpolateCrust ( Mesh &msh, Discontinuity &dis );
  void interpolateTopo  ( Mesh &msh, Discontinuity &dis );
  void interpolate      ( Mesh &msh, Model_file &mod, Discontinuity &dis );  
  int  recover      ( double &testX, double &testY, double &testZ, Mesh &msh,
                      double &c11, double &c12, double &c13, double &c14, 
                      double &c15, double &c16, double &c22, double &c23, 
                      double &c24, double &c25, double &c26, double &c33, 
                      double &c34, double &c35, double &c36, double &c44, 
                      double &c45, double &c46, double &c55, double &c56, 
                      double &c66, double &rho, char mode ); 
                      
private:
  double taper      ( double &x, double &y, double &z, Model_file &mod );
  
};

class Attenuation
{

public:
  
  int nRelaxationMechanisms = 3;
  
  double freqRef = 1;
  double *tau_s  = new double [3];
  double *D      = new double [3];
  
  std::string qModelName = "QL6";
  
  double QL6      ( double &rad );
  double correct ( std::string &model, double &rad );


};