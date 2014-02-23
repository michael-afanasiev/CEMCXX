#include "classes.hpp"
#include <time.h>
#include <cmath>
#include "kdtree.h"

using namespace std;

int main () 
{
  
  int i;
  
  ifstream myfile;
  
  Region        reg;
  
  Constants     con;
  Exodus_file   exo;
  Driver        drv;
  Model_file    mod;
  Utilities     util;
  Discontinuity dis;
  Interpolator inter;
  
  cout << "Begin model building.\n";  
  
  // ********************************************************************* //
  //                       READING THE PARAMETER FILE                      //
  // ********************************************************************* //
    
  // Determine how many lines to read from parameter file.
  drv.openDriver ( myfile );  
  
  // ********************************************************************* //
  //                            READ IN PARAMS                             //
  // ********************************************************************* //
  
  mod.populateParams ( drv, exo );        

  // ********************************************************************* //
  //                            READ IN MODEL                              //
  // ********************************************************************* //

  
  mod.read ();  
  dis.read ();
  
  if ( mod.rotAng != 0 ) {
    util.rotate ( mod );
  }  
  
  mod.openUp          ();
  
  
  // ********************************************************************* //
  //                           CONSTRUCT THE MESH                          //
  // ********************************************************************* //
  
  exo.merge                ( reg, mod );
  
  if ( mod.intentions == "INTERPOLATE" ) {
    
    cout << "\n----- Interpolating -----\n";
    
    mod.createKDTreeUnpacked ( );
    dis.createKDTreePacked   ( );
    
    for ( vector < Exodus_file > :: iterator it=reg.regionsExo.begin();
      it!=reg.regionsExo.end(); ++it ) {  
        
      cout << "\n";

      Mesh msh;
      Interpolator ipl;    
    
      it -> openFile     ( it -> fname );        
      msh.getInfo        ( it -> idexo, 'p' );    
      ipl.interpolate    ( msh, mod, dis );
      it -> writeParams  ( msh );
      msh.deallocateMesh ();
      it -> closeFile    ();      
    
    }
  }
  
  if ( mod.intentions == "EXTRACT" ) {
    
    cout << "\n----- Extracting -----\n";
        
    for ( vector < Exodus_file > :: iterator it=reg.regionsExo.begin();
      it!=reg.regionsExo.end(); ++it ) {
        
      cout << "\n";
        
      double c11, c12, c13, c14, c15, c16,c22, c23, c24, c25, c26, c33, c34;
      double c35, c36, c44, c45, c46, c55, c56, c66, rho, testX, testY, testZ;
    
      Mesh         msh;
      Interpolator ipl;
    
      it -> openFile           ( it -> fname );
      msh.getInfo              ( it -> idexo, 'p' );
      msh.getConnectivity      ( it -> idexo );
      msh.createKDTreeUnpacked ( );            
      
      int eSum = 0;
      for ( int r=0; r<mod.col_deg.size(); r++ ) {
        for ( int i=0; i<mod.col_rad[r].size(); i++ ) {
          for ( int j=0; j<mod.lon_rad[r].size(); j++ ) {
            for ( int k=0; k<(mod.rad[r].size()-1); k++ ) {
                            
              int ind = (i * mod.lon_rad[r].size() * (mod.rad[r].size()-1) )
                + (j * (mod.rad[r].size()-1) ) + ( k ) + eSum;
                            
              if ( (mod.rad[r][k]     <= msh.radMax) && 
                   (mod.rad[r][k]     >= msh.radMin) &&
                   (mod.lon_rad[r][j] <= msh.lonMax) &&
                   (mod.lon_rad[r][j] >= msh.lonMin) &&
                   (mod.col_rad[r][i] <= msh.colMax) &&
                   (mod.col_rad[r][i] >= msh.colMin ) ) 
              {
                            
                util.colLonRadRad2xyz ( mod.col_rad[r][i], mod.lon_rad[r][j], 
                  mod.rad[r][k], testX, testY, testZ );
                                          
                int pass = ipl.recover ( testX, testY, testZ, msh.tree, msh,
                  c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26,
                  c33, c34, c35, c36, c44, c45, c46, c55, c56, c66,
                  rho, 'p' );                  
                  
                if ( mod.input_model_physics == "TTI" ) {
                  mod.vsh[r][ind] = sqrt (c44 / rho);
                  mod.vsv[r][ind] = sqrt (c55 / rho);
                  mod.vpp[r][ind] = sqrt (c22 / rho);
                  mod.rho[r][ind] = rho;
                }                
              }
            }
          }
          
          cout << "                                                           " 
            << "\r" << std::flush;
          cout << "CoLat loops left for this chunk & region: " << 
            mod.col_rad[r].size() - (i + 1) << "\r" << std::flush;                    
        }
        
        cout << "\n";
        eSum += mod.col_rad[r].size() * mod.lon_rad[r].size() * 
          ( mod.rad[r].size() - 1 );
      }
      
      msh.deallocateMesh ();
      it -> closeFile    ();      
    }
    
    mod.writeSES3D ();           
  }
                  
  return 0;
    
}
