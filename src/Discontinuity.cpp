#include <sstream>
#include <cmath>

#include "kdtree.h"
#include "classes.hpp"
void Discontinuity::read ()
{
  
  Model_file mod;
  
  std::string cmd = "./dat/discontinuities";
  int dum1;

  // Read moho depth and crustal parameters.
  mod.populateSES3D ( cmd + "/crust_x_smooth",   dum1, crust_col_deg, 
  'c' );
  mod.populateSES3D ( cmd + "/crust_y_smooth",   dum1, crust_lon_deg, 
  'c' );
  mod.populateSES3D ( cmd + "/crust_vs_smooth",  dum1, crust_vs, 'p' );
  mod.populateSES3D ( cmd + "/crust_dep_smooth", dum1, crust_dp, 'p' );  

  for ( size_t i=0; i<crust_lon_deg[0].size(); i++ ) {
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
  crust_col_deg_unpack.reserve ( crust_col_rad[0].size()*
    crust_lon_rad[0].size() );
  crust_lon_deg_unpack.reserve ( crust_col_rad[0].size()*
    crust_lon_rad[0].size() );
  
      
    
  for ( size_t r=0; r!=crust_col_rad.size(); r++ ) {
    for ( size_t i=0; i!=crust_col_rad[r].size(); i++ ) {
      for ( size_t j=0; j!=crust_lon_rad[r].size(); j++ ) {
    
        KDdatCrust[l] = l;
        kd_insert3 ( crustTree, crust_col_deg[r][i], crust_lon_deg[r][j], 
          con.R_EARTH, &KDdatCrust[l] );
        crust_col_deg_unpack[l] = crust_col_deg[r][i];
        crust_lon_deg_unpack[l] = crust_lon_deg[r][j];
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
  for ( size_t i=0; i<lonElv.size(); i++ ) 
  { 
    KDdatElv[i] = i;
    {   
      kd_insert3 ( elvTree, colElv[i], lonElv[i], con.R_EARTH, &KDdatElv[i] );
    }
  }
  
}

void Discontinuity::lookCrust ( Mesh &msh, double &mshCol, double &mshLon, 
                                double &mshRad, int &mshInd, bool &checkCrust,
                                bool &smoothCrust, double &upTap, 
                                double &downTap, Model_file &mod )
{
  
  Constants con;
  double rho, vpv;
  
  if ( mshRad > (con.R_EARTH - 100) ) 
  {          
    kdres *set = kd_nearest3      ( crustTree, mshCol, mshLon, con.R_EARTH );         
    void *ind  = kd_res_item_data ( set );
    int  point = * ( int * ) ind;
    
    kd_res_free (set);
    
    /* The moho is defined in a weird way ( depth from sea level if in the 
    ocean, and depth from elevation if in the crust). First, convert crust 
    elevation to km, and then decided whether we're taking the sea level or
    or crustial surface as reference */
    double ref;
    if ( msh.elv[mshInd] <= 0. )
    {
      ref = con.R_EARTH;
    }
    else
    {
      ref = con.R_EARTH + msh.elv[mshInd] / 1000.;
    }
    
    // Do a bilinear interpolation on the Depth.
    double interpDep;
    double interpVs;
    getCrustDepth ( mshCol, mshLon, point, interpDep, "dep" );
    getCrustDepth ( mshCol, mshLon, point, interpVs,  "vel" );
    

    /* Here get crustal thickness for tapering purposes (differs if we're in)
    oceanic of continental crust. This is the thickness of the crust, minus an 
    ocean layer. */
    double crustThick = 0.;
    if ( msh.elv[mshInd] <= 0. )
    {
      crustThick = interpDep + ( msh.elv[mshInd] );
    }
    else
    {
      crustThick = interpDep;
    }
    
    /* Figure out whether we need to smooth the crust */
    smoothCrust = false;
    smoothCosine ( crustThick, interpDep, mshRad, downTap, upTap, smoothCrust );
    
    /* If we're in the crust */
    // This used to say if mshRad >= (ref - interpDep) as well.
    if ( mod.intentions == "CRUST" ) 
    {
      double crust_vsv = interpVs - con.aniCorrection;
      double crust_vsh = interpVs;
      
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
      msh.c11[mshInd] = upTap * A   * downTap * msh.c11[mshInd];
      msh.c22[mshInd] = upTap * A   * downTap * msh.c22[mshInd];
      msh.c33[mshInd] = upTap * A   * downTap * msh.c33[mshInd];        
      msh.c12[mshInd] = upTap * F   * downTap * msh.c12[mshInd];
      msh.c13[mshInd] = upTap * F   * downTap * msh.c13[mshInd];
      msh.c23[mshInd] = upTap * S   * downTap * msh.c23[mshInd];
      msh.c44[mshInd] = upTap * N   * downTap * msh.c44[mshInd];
      msh.c55[mshInd] = upTap * L   * downTap * msh.c55[mshInd];
      msh.c66[mshInd] = upTap * L   * downTap * msh.c66[mshInd];      
      msh.rho[mshInd] = upTap;// * rho * downTap * msh.rho[mshInd];        
    }
    
    if ( mshRad >= (ref - interpDep) )
    {
      checkCrust = true;
    }
  }    
  
}

void Discontinuity::getCrustDepth ( double &mshCol, double &mshLon, int &point,
                                    double &par, std::string mode )
{
  
  Constants con;
  
  // Stride of one degree in crust07.
  int degStri=2;
  
  // Initializations.
  int lonDir = 0; 
  int colDir = 0;
  
  // Create vector to hold 4 nodes of smoothing square.
  std::vector <int> nodes;
  nodes.reserve (4);
    
  // Determine on which side of the requested col. the closest point is.
  if ( mshCol > crust_col_deg_unpack[point] )
    colDir = 1;  
  if ( mshCol <= crust_col_deg_unpack[point] )
    colDir = -1;
  
  // Determine on which side of the requested lon. the closest point is.  
  if ( mshLon > crust_lon_deg_unpack[point] )
    lonDir = 1;
  if ( mshLon <= crust_lon_deg_unpack[point] )
    lonDir = -1;

  // Create the KDTree requests for four points of the interpolating square.
  double col1 = crust_col_deg_unpack[point];
  double col2 = crust_col_deg_unpack[point] + colDir * degStri;
  double lon1 = crust_lon_deg_unpack[point];
  double lon2 = crust_lon_deg_unpack[point] + lonDir * degStri;
  
  // Fix col. wrapping.
  if ( col2 > 180 )
    col2 = 180 - ( col1 - 180 );
  if ( col2 < 0 )
    col2 = 0;
  
  // Fix lon. wrapping.
  if ( lon2 > 180 )
    lon2 = lon2 - 360.;
  if ( lon2 < -180 )
    lon2 = lon2 + 360.;    
  
  // Find the four points of the rectangle in spherical space.
  kdres *set1 = kd_nearest3 ( crustTree, col1, lon1, con.R_EARTH );
  kdres *set2 = kd_nearest3 ( crustTree, col1, lon2, con.R_EARTH );
  kdres *set3 = kd_nearest3 ( crustTree, col2, lon1, con.R_EARTH );
  kdres *set4 = kd_nearest3 ( crustTree, col2, lon2, con.R_EARTH );
  
  // Extract the indices from the result structures.
  void *ind1 = kd_res_item_data ( set1 ); 
  void *ind2 = kd_res_item_data ( set2 );
  void *ind3 = kd_res_item_data ( set3 );
  void *ind4 = kd_res_item_data ( set4 );
  
  // Build the square based on the four possible quadrants.
  if ( colDir == 1 && lonDir == (-1) )
  {
    nodes[0] = * ( int * ) ind2;
    nodes[1] = * ( int * ) ind1;
    nodes[2] = * ( int * ) ind3;
    nodes[3] = * ( int * ) ind4;
  } 
  else if ( colDir == 1 && lonDir == 1 )
  {
    nodes[0] = * ( int * ) ind1;
    nodes[1] = * ( int * ) ind2;
    nodes[2] = * ( int * ) ind4;
    nodes[3] = * ( int * ) ind3;  
  }
  else if ( colDir == (-1) && lonDir == 1 )
  {
    nodes[0] = * ( int * ) ind3;
    nodes[1] = * ( int * ) ind4;
    nodes[2] = * ( int * ) ind2;
    nodes[3] = * ( int * ) ind1;
  }
  else if ( colDir == (-1) && lonDir == (-1) )
  {
    nodes[0] = * ( int * ) ind4;
    nodes[1] = * ( int * ) ind3;
    nodes[2] = * ( int * ) ind1;
    nodes[3] = * ( int * ) ind2;
  }            
  
  // Deallocate memory for results.
  kd_res_free ( set1 );
  kd_res_free ( set2 );
  kd_res_free ( set3 );
  kd_res_free ( set4 );
  
  // Conditionally build the weights depending on the quandrants selected above.
  double t = 0;
  double u = 0;
  if ( colDir == 1 )
    t = ( mshCol - col1 ) / ( col2 - col1 );
  if ( colDir == (-1) )                   
    t = ( mshCol - col2 ) / ( col1 - col2 );
  if ( lonDir == 1 )
    u = ( mshLon - lon1 ) / ( lon2 - lon1 );
  if ( lonDir == (-1) )
    u = ( mshLon - lon2 ) / ( lon1 - lon2 );
    
  /* Get the bilinearly interpolated value. Inspired by Numerical recipes,
  but with a swwitched axis. */
  if ( mode == "dep" )
  {
    par = (1 - t) * (1 - u) * crust_dp[0][nodes[0]] +
      t * (1 - u) * crust_dp[0][nodes[3]] +
      t * u * crust_dp[0][nodes[2]] +
      (1 - t) * u * crust_dp[0][nodes[1]];  
  }
    
  if ( mode == "vel" )
  {
    par = (1 - t) * (1 - u) * crust_vs[0][nodes[0]] +
      t * (1 - u) * crust_vs[0][nodes[3]] +
      t * u * crust_vs[0][nodes[2]] +
      (1 - t) * u * crust_vs[0][nodes[1]];  
  }
  

#ifdef VISUAL_DEBUG
  std::ofstream myfile;
  myfile.open ( "surfInterp.txt", std::ios::out );
  
  myfile << crust_col_deg_unpack[nodes[0]] << " " 
    << crust_lon_deg_unpack[nodes[0]] << " " << crust_dp[0][nodes[0]] 
    << std::endl;
      
  myfile << crust_col_deg_unpack[nodes[1]] << " " 
    << crust_lon_deg_unpack[nodes[1]] << " " << crust_dp[0][nodes[1]] 
    << std::endl;
        
  myfile << crust_col_deg_unpack[nodes[2]] << " " 
    << crust_lon_deg_unpack[nodes[2]] << " " << crust_dp[0][nodes[2]] 
    << std::endl;
      
  myfile << crust_col_deg_unpack[nodes[3]] << " " 
    << crust_lon_deg_unpack[nodes[3]] << " " << crust_dp[0][nodes[3]] 
    << std::endl;
  
  myfile << mshCol << " " 
    << mshLon << " " << dep
    << std::endl;
  
  std::cout << colDir << " " << lonDir << std::endl;
  std::cout << mshCol << " " << col1 << " " << col2 << " " << std::endl;
  std::cout << mshLon << " " << lon1 << " " << lon2 << " " << std::endl;  
  std::cout << t << " " << u << std::endl;
  std::cout << crust_dp[0][nodes[0]] << " " << crust_dp[0][nodes[1]] << " " 
    << crust_dp[0][nodes[2]] << " " << crust_dp[0][nodes[3]] << std::endl;
  std::cout << dep << std::endl;
  
  myfile.close();
  std::cin.get();     
  
  std::cout << colDir << " " << lonDir << std::endl;
  std::cout << mshCol << " " << col1 << " " << col2 << " " << std::endl;
  std::cout << mshLon << " " << lon1 << " " << lon2 << " " << std::endl;  
  std::cout << t << " " << u << std::endl;
  std::cout << crust_dp[0][nodes[0]] << " " << crust_dp[0][nodes[1]] << " " 
    << crust_dp[0][nodes[2]] << " " << crust_dp[0][nodes[3]] << std::endl;
  std::cout << dep << std::endl;
#endif
  
}

void Discontinuity::smoothCosine ( double &crustThick, double &centerDep, 
                                   double &rad, double &downTap, 
                                   double &upTap, bool &smoothTrue )
{
  
  Constants con;
  
  smoothTrue           = false;
  double maxSmooth     = 15.;
  double minSmooth     = 3.;
  double smoothPercent = 0.15;

  /* Smooth at most over smoothPercent of the crust. If we're below the min
  or above the max, don't smooth any more/less */
  double smoothDist = smoothPercent * crustThick;
  if ( smoothDist < minSmooth )
    smoothDist = minSmooth;
  if ( smoothDist > maxSmooth )
    smoothDist = maxSmooth;
  
  /* If we're in the smoothing region, apply a cosine smoother. UpTap slowly moves
   * to crustal values */
  double disconRad = con.R_EARTH - centerDep;
  if ( rad <= disconRad + (smoothDist / 2.) &&
       rad >= disconRad - (smoothDist / 2.) )
  {
    double length  = rad - (disconRad - (smoothDist / 2.));
    downTap = (1/2.) * (cos ( con.PI * length / smoothDist ) + 1);
    upTap   = 1 - downTap;
    
    smoothTrue = true;
    
    if ( upTap < 0. || upTap > 1. )
    {
      std::cout << "Error in discontinuity smoothing" << std::endl;
      exit (EXIT_FAILURE);
    }
    if ( downTap < 0. || downTap > 1. )
    {
      std::cout << "Error in discontinuity smoothing" << std::endl;
      exit (EXIT_FAILURE);
    }
          
  }

  /* If we're above the smoothing region, it's all crust baby */
  if ( rad > disconRad + (smoothDist / 2.) )
  {
    upTap   = 1.;
    downTap = 0.;
  }

  /* And if we're below, it's all S20 */
  if ( rad < disconRad - (smoothDist / 2.) )
  {
    upTap   = 0.;
    downTap = 1.;
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
