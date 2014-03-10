#include "kdtree.h"
#include "classes.hpp"

void Discontinuity::read ()
{
  
  Model_file mod;
  
  std::string cmd = "./dat/discontinuities";
  int dum1, dum2;

  // Read moho depth and crustal parameters.
  mod.populateSES3D ( cmd + "/crust_x_smooth",   dum1, crust_col_deg, 
  'c' );
  mod.populateSES3D ( cmd + "/crust_y_smooth",   dum1, crust_lon_deg, 
  'c' );
  mod.populateSES3D ( cmd + "/crust_vs_smooth",  dum1, crust_vs, 'p' );
  mod.populateSES3D ( cmd + "/crust_dep_smooth", dum1, crust_dp, 'p' );  

  for ( int i=0; i<crust_lon_deg[0].size(); i++ ) {
    if ( crust_lon_deg[0][i] > 180. )
      crust_lon_deg[0][i] = crust_lon_deg[0][i] - 360.;      
  }

  mod.populateRadians ( crust_col_deg, crust_col_rad );
  mod.populateRadians ( crust_lon_deg, crust_lon_rad );
  

}

void Discontinuity::createKDTreePacked () 
{
  
  Constants con;
    
  std::cout << "Creating KDTree ( crust ).\n";
  
  tree = kd_create (3);

  int l    = 0;
  int *dat = new int [crust_col_rad[0].size()*crust_lon_rad[0].size()]();
    
  for ( int r=0; r!=crust_col_rad.size(); r++ ) {
    for ( int i=0; i!=crust_col_rad[r].size(); i++ ) {
      for ( int j=0; j!=crust_lon_rad[r].size(); j++ ) {
    
        dat[l] = l;
        kd_insert3 ( tree, crust_col_deg[r][i], crust_lon_deg[r][j], 
          con.R_EARTH, &dat[l] );
        l++;
    
      }
    }               
  }
  
}

void Discontinuity::lookCrust ( Mesh &msh, double &mshCol, double &mshLon, 
                                double &mshRad, int &mshInd )
{
  
  Constants con;
  
  if ( mshRad > (con.R_EARTH - 100) ) {          

    kdres *set = kd_nearest3      ( tree, mshCol, mshLon, con.R_EARTH );         
    void *ind  = kd_res_item_data ( set );
    int  point = * ( int * ) ind;
    
    if ( mshRad >= (con.R_EARTH - crust_dp[0][point]) ) {
  
      double crust_vsv = crust_vs[0][point] - con.aniCorrection;
      double crust_vsh = crust_vs[0][point];
        
      double N = msh.rho[mshInd] * crust_vsh * crust_vsh;
      double L = msh.rho[mshInd] * crust_vsv * crust_vsv;
        
      double A = msh.c22[mshInd];
      double S = A - 2 * N;
      double F = A - 2 * L;
  
      inCrust         = true;
      msh.c12[mshInd] = F;
      msh.c13[mshInd] = F;
      msh.c23[mshInd] = S;
      msh.c44[mshInd] = N;
      msh.c55[mshInd] = L;
      msh.c66[mshInd] = L;                
  
    }
  }    
  
}
