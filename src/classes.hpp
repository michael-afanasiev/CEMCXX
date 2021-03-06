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
  std::vector < std::string > colReg;
  std::vector < std::string > lonReg;


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
  double tiny          = 1;
  double bigtiny       = 3;
  double tinytiny      = 1e-2;
  double oneDegRad     = 1. * PI / o80;
  
  double innerCoreRad = 1221.0;
  double outerCoreRad = 3480.0;
  double RTHO         = 5371.0;
  double R670         = 5701.0;
  double R400         = 5971.0; 
  double R100         = 6271.0;
  double R052         = 6319.0;
  double R020         = 6351.0;
  
};

class Model_file
{
  
public:
  
  void readSES3D     ();
  void readSPECFEM3D ();
  
  void populateSES3D ( std::string name, int &num_regions, 
    std::vector <std::vector <double> > &vec , char ftype );
      
  int *region;
  int *KDdat1;
  int *KDdat2;
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
  std::vector <double> lonUnwrap;
  std::vector <double> colUnwrap;
  
  std::vector <short> r;
  
  std::string specFileName;
      
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
  bool radReg1    = false;
  bool radReg2    = false;
  bool radReg3    = false;
  bool radReg4    = false;
  bool radReg5    = false;
  bool radReg6    = false;
  bool radReg7    = false;
  bool radReg8    = false;
  bool radReg9    = false;
  bool radReg10   = false;
  bool radReg11   = false;
  bool radReg12   = false;
  bool radReg13   = false;
  bool radReg14   = false;
  bool radReg15   = false;
  bool radReg16   = false;
  bool radReg17   = false;
  bool radReg18   = false;
  bool radReg19   = false;
  bool radReg20   = false;
  bool radReg21   = false;
  bool radReg22   = false;
  bool radReg23   = false;
  bool radReg24   = false;
  bool radReg25   = false;
  bool radReg26   = false;
  bool radReg27   = false;
  bool radReg28   = false;
  bool radReg29   = false;
  bool radReg30   = false;
  bool radReg31   = false;
  bool radReg32   = false;
  bool radReg33   = false;
  bool radReg34   = false;
  bool radReg35   = false;
  bool radReg36   = false;
  bool radReg37   = false;
  bool radReg38   = false;
  bool radReg39   = false;
  bool radReg40   = false;
  bool radReg41   = false;
  bool radReg42   = false;
  bool radReg43   = false;
  bool radReg44   = false;
  bool radReg45   = false;
  bool wrapAround = false;
  bool doRotate   = false;
  bool kernel1d   = false;
  bool kernel3d   = false;
  bool overwriteCrust;

  
  kdtree *tree1;
  kdtree *tree2;

  std::vector<kdtree*> treeVec;
  std::vector <int> kdRegions;

  std::vector <std::vector <double> > col_rad;
  std::vector <std::vector <double> > lon_rad;
  std::vector <std::vector <double> > col_deg;
  std::vector <std::vector <double> > lon_deg;
  std::vector <std::vector <double> > rad;
  std::vector <std::vector <double> > vsh;
  std::vector <std::vector <double> > vsv;
  std::vector <std::vector <double> > rho;
  std::vector <std::vector <double> > vpp;  
  
  std::vector <double> minRadReg;
  std::vector <double> maxRadReg;

  std::string input_model_directory;
  std::string input_model_file_type;
  std::string input_model_physics;
  std::string intentions;
  std::string output_model_physics;
  std::string mesh_directory;
  std::string subset;
  std::string terraFileName;
  std::string terraFileProc;

  // Internal functions.  
  void read                   ( );
  void colLonRad2xyzSES3D     ( );
  void getMinMaxRegionSES3D   ( );
  void populateRadiansSES3D   ( );
  void openUp                 ( );
  void projectSubspaceSES3D   ( );
  void writeSES3D             ( );
  void deallocate             ( );
  void readTERRAGRID          ( );
  void writeTERRAGRID         ( );
  void projectSubspaceSPECFEM ( );
  void getSpecFileName        ( int &, int & );
  void getTerraFileName       ( int & );
  void createKDTreeUnpacked   ( Mesh &msh );
  int  writeNetCDF            ( std::vector <std::vector <double>> &par,
                                std::string name );
  void populateParams         ( Driver &drv, Exodus_file &exo );
  void dePopulateSES3D        ( std::string, std::vector<std::vector<double>> );
  void populateRadians        ( std::vector < std::vector <double> > &deg, 
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
  
  std::vector <double> crust_col_deg_unpack;
  std::vector <double> crust_lon_deg_unpack;
  
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
  void smoothCosine         ( double &crustThick, double &centerDep, 
                              double &rad, double &downTap, double &upTap, 
                              bool &smoothTrue );
  void lookCrust            ( Mesh &msh, double &mshCol, double &mshLon, 
                            double &mshRad, int &mshInd, bool &checkCrust,
                            bool &smoothCrust, double &upTap, double &downTap, Model_file &mod );
  void lookTopo             ( Mesh &msh, double &mshCol, double &mshLon, 
                              double &mshRad, int &mshInd );
                              
private:
  void getCrustDepth ( double &, double &, int &, double &, std::string );

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
  void checkRegionExtr  ( double x, double y, double z, short r,
                          double &xUse, double &yUse, double &zUse );
  void checkMeshEdge    ( double &col, double &lon, Mesh &msh ); 
  void pullInRad        ( double &col, double &lon, double &rad,
                          double &x,   double &y,   double &z,
                          Mesh &msh );
  int getFilesize       ( std::string fname );

  void pullRad          ( double &col, double &lon, double &rad, Mesh &msh, bool &fullSearch );
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
  
  void readDriver     ( std::ifstream &myfile );
  void checkUsage     ( Model_file &mod, std::string mode );
  void initialize     ( Model_file &mod, Discontinuity &dis, Utilities &utl, 
                        Exodus_file &exo, Region &reg );
  void populateParams ( Model_file &mod );
  void report         ( Model_file &mod );
  
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
  int *sideSet;
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
  
  bool colReg000_090 = false;
  bool colReg090_180 = false;
  bool lonReg000_090 = false;
  bool lonReg090_180 = false;
  bool lonReg180_270 = false;
  bool lonReg270_360 = false;

  kdtree *tree;
  
  std::multimap <int, std::vector <int> > elemOrder;  
  std::vector < std::vector <int> > refineElemConn;

  std::vector <int> masterElemConn;  
  char name;
  char title [MAX_LINE_LENGTH+1]; 
  char elem_type [MAX_LINE_LENGTH+1];
  
  std::string regString;
  
  // Internal functions.    
  void populateCoord          ( );
  void populateParams         ( );
  void allocateMesh           ( );
  void createKDTreeUnpacked   ( );
  void getMinMaxRad           ( );
  void getInfo                ( int exoid );
  void getSideSet             ( int exoid );
  void getConnectivity        ( int exoid );
  void getNodeNumMap          ( int exoid );
  void getElemNumMap          ( int exoid );
  void deallocateMesh         ( Model_file &mod );
  void getRegion              ( Region &reg, int fileIter );
           
};

class Interpolator
{
  
public:
 
  int fallBackCount = 10000;

  std::vector <int> refineArr;
  std::vector < std::vector <int> > elemWithin;
    
  void interpolateCrust ( Mesh &msh, Discontinuity &dis, Model_file &mod );
  void interpolateTopo  ( Mesh &msh, Discontinuity &dis );
  void interpolate      ( Mesh &msh, Model_file &mod, Discontinuity &dis );  
  int  recover          ( double &testX, double &testY, double &testZ, 
                          Mesh &msh,
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
