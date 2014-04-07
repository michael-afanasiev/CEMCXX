#include <sstream>
#include <cmath>

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
  
  crustTree = kd_create (3);

  int l    = 0;
  KDdatCrust = new int [crust_col_rad[0].size()*crust_lon_rad[0].size()]();
    
  for ( int r=0; r!=crust_col_rad.size(); r++ ) {
    for ( int i=0; i!=crust_col_rad[r].size(); i++ ) {
      for ( int j=0; j!=crust_lon_rad[r].size(); j++ ) {
    
        KDdatCrust[l] = l;
        kd_insert3 ( crustTree, crust_col_deg[r][i], crust_lon_deg[r][j], 
          con.R_EARTH, &KDdatCrust[l] );
        l++;
    
      }
    }               
  }
  
}

void Discontinuity::createKDTreeUnpacked ( )
{
  
  Constants con;
  
  std::cout << "Creating KDTree ( model ).\n";
  elvTree  = kd_create (3);
  KDdatElv = new int [ lonElv.size() ];

#pragma omp parallel for
  for ( int i=0; i<lonElv.size(); i++ ) 
  { 
    KDdatElv[i] = i;
    {   
      kd_insert3 ( elvTree, colElv[i], lonElv[i], con.R_EARTH, &KDdatElv[i] );
    }
  }
  
}

void Discontinuity::lookCrust ( Mesh &msh, double &mshCol, double &mshLon, 
                                double &mshRad, int &mshInd, bool &checkCrust )
{
  
  Constants con;
  double rho, vpv;
  
  if ( mshRad > (con.R_EARTH - 100) ) {          

    kdres *set = kd_nearest3      ( crustTree, mshCol, mshLon, con.R_EARTH );         
    void *ind  = kd_res_item_data ( set );
    int  point = * ( int * ) ind;
    
    /* The moho is defined in a weird way ( depth from sea level if in the 
    ocean, and depth from elevation if in the crust). First, convert crust 
    elevation to km, and then decided whether we're taking the sea level or
    or crust as reference */
    double ref;
    if ( msh.elv[mshInd] <= 0. )
    {
      ref = con.R_EARTH;
    }
    else
    {
      ref = con.R_EARTH + msh.elv[mshInd] / 1000.;
    }
    
    if ( mshRad >= (ref - crust_dp[0][point]) ) 
    {
      double crust_vsv = crust_vs[0][point] - con.aniCorrection;
      double crust_vsh = crust_vs[0][point];
      
      // Scaling from isotropic vs to rho (Fichtner, multiscale)
      rho = 0.2277 * crust_vsh + 2.016;
      
      // Scaling from vs to vp
      vpv = 1.5399 * crust_vsh + 0.840;        
        
      double N = rho * crust_vsh * crust_vsh;
      double L = rho * crust_vsv * crust_vsv;
        
      double A = rho * vpv * vpv;
      double S = A - 2 * N;
      double F = A - 2 * L;
  
      checkCrust      = true;
      msh.c11[mshInd] = A;
      msh.c22[mshInd] = A;
      msh.c33[mshInd] = A;        
      msh.c12[mshInd] = F;
      msh.c13[mshInd] = F;
      msh.c23[mshInd] = S;
      msh.c44[mshInd] = N;
      msh.c55[mshInd] = L;
      msh.c66[mshInd] = L;                
      msh.rho[mshInd] = rho;  
    }
  }    
  
}

void Discontinuity::lookTopo ( Mesh &msh, double &mshCol, double &mshLon, 
                               double &mshRad, int &mshInd )
{
  
  Constants con;
  
  if ( mshRad > ( con.R_EARTH - 100 ) )
  {
    kdres *set = kd_nearest3      ( elvTree, mshCol, mshLon, con.R_EARTH );         
    void *ind  = kd_res_item_data ( set );
    int  point = * ( int * ) ind;
    
    msh.elv[mshInd] = elv[point];
  }
  
}

void Discontinuity::readTopography ( )
{
  
  std::ifstream myfile;
  Constants con;
  
  std::cout << "Reading topography." << std::endl;
  myfile.open ( "./dat/discontinuities/10MinuteTopoGrid.txt", std::ios::in );
  
  std::string line;
  while ( std::getline (myfile, line) )    
  {
    
    std::stringstream linestream (line);
    std::string value;
    
    int i = 0;
    while ( std::getline (linestream, value, ',') )
    {
      if ( i == 0 )
      {
        double lon = stod ( value );
        if ( lon < 0. )
          lon = 360 + lon;
        lonElv.push_back ( stod (value) );
      }        
      if ( i == 1 )
        colElv.push_back ( 90. - stod (value) );
      if ( i == 2 )
        elv.push_back ( stod (value) );
      i++;
    }     
  }
        
}

void Discontinuity::deallocate ( )
{
  
  kd_free ( crustTree );
  kd_free ( elvTree );
  delete [] KDdatCrust;
  delete [] KDdatElv;
  
}
