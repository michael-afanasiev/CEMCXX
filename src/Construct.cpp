#include "classes.hpp"
#include <time.h>
#include <cmath>
#include "kdtree.h"

using namespace std;

int main () 
{
  
  ifstream myfile;
  
  Region        reg;  
  Constants     con;
  Exodus_file   exo;
  Driver        drv;
  Model_file    mod;
  Utilities     utl;
  Discontinuity dis;
  
  cout << "Begin model building.\n";  
      
  // Read parameter file.
  drv.openDriver     ( myfile );  
  
  // Copy parameters to model file.  
  mod.populateParams ( drv, exo );        
  
  // Read model file.
  mod.read           ( );  
  
  // Read discontinuities.
  dis.read           ( );
  
  // Rotate model coordinates if necessary.
  utl.inquireRotate  ( mod );
      
  // Project to elastic tensor.
  mod.openUp         ( );
  
  // Determine which exodus files to read.  
  exo.merge          ( reg, mod );
  
  if ( mod.intentions == "INTERPOLATE" )   
  {
    
    cout << "\n----- Interpolating -----\n";
    cout << "\n";
    
    mod.createKDTreeUnpacked ( );
    dis.createKDTreePacked   ( );
    
    for ( vector < Exodus_file > :: iterator exoFile=reg.regionsExo.begin();
      exoFile!=reg.regionsExo.end(); ++exoFile )       
    {  
        
      cout << "\n";

      Mesh msh;
      Interpolator ipl;    
    
      exoFile -> openFile    ( exoFile -> fname );        
      msh.getInfo            ( exoFile -> idexo, 'p' );    
      ipl.interpolate        ( msh, mod, dis );
      exoFile -> writeParams ( msh );
      msh.deallocateMesh     ( mod );
      exoFile -> closeFile   ( );      
    
    }
  }
  
  if ( mod.intentions == "EXTRACT" ) 
  {
    
    cout << "\n----- Extracting -----\n";
        
    for ( vector < Exodus_file > :: iterator exoFile=reg.regionsExo.begin();
      exoFile!=reg.regionsExo.end(); ++exoFile ) 
    {
        
      cout << "\n";
        
      double c11, c12, c13, c14, c15, c16,c22, c23, c24, c25, c26, c33, c34;
      double c35, c36, c44, c45, c46, c55, c56, c66, rho, testX, testY, testZ;
      double col, lon, rad;
    
      Mesh         msh;
      Interpolator ipl;
    
      exoFile -> openFile      ( exoFile -> fname );
      msh.getInfo              ( exoFile -> idexo, 'p' );
      msh.getConnectivity      ( exoFile -> idexo );
      msh.createKDTreeUnpacked ( );            
      
      cout << "Extracting." << endl;
      for ( int i=0; i<mod.x.size(); i++ )         
      {
        
        
        utl.rotateForward    ( mod.x[i], mod.y[i], mod.z[i], testX, 
                               testY, testZ, mod );
        
        utl.xyz2ColLonRadRad ( testX, testY, testZ, col, lon, rad ); 
                                      
                              
        if ( (rad <= msh.radMax) && 
             (rad >= msh.radMin) &&
             (lon <= msh.lonMax) &&
             (lon >= msh.lonMin) &&
             (col <= msh.colMax) &&
             (col >= msh.colMin ) ) 
        {
                                               
                                               
        int pass = ipl.recover ( testX, testY, testZ, msh,  c11, c12, c13, c14, 
          c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, 
          c55, c56, c66, rho, 'p' ); 
        
        mod.c11[i]    = c11;
        mod.c12[i]    = c12;
        mod.c13[i]    = c13;
        mod.c22[i]    = c22;
        mod.c23[i]    = c23;
        mod.c33[i]    = c33;
        mod.c44[i]    = c44;
        mod.c55[i]    = c55;
        mod.c66[i]    = c66;            
        mod.rhoUnwrap[i] = rho;   
        
        }                          
      }        

      msh.deallocateMesh   ( mod );
      exoFile -> closeFile ( );                 
    }
    
    mod.projectSubspace ( );
    mod.writeSES3D      ( );           
  }
                  
  return 0;    
}