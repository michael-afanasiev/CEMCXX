#include <iostream>
#include <cmath>
#include <assert.h>
#include <ctype.h>
#include <time.h>

#include "kdtree.h"
#include "classes.hpp"

using namespace std;

void Interpolator::interpolate ( Mesh &msh, Model_file &mod ) 
{

  Utilities util;
  
  // Create KDTree.
  cout << "Creating KDTree.\n";
  kdtree *tree = kd_create (3);  
  int *dat     = new int [mod.num_p];
  for ( int i=0; i<mod.num_p; i++ ) {
    dat[i] = i;
    kd_insert3 ( tree, mod.x[i], mod.y[i], mod.z[i], &dat[i] );
  }      
  
  cout << "Interpolating.\n";
  
  for ( int i=0; i<msh.num_nodes; i++ ) {
    
    double mshCol;
    double mshLon;
    double mshRad;
    
    util.xyz2ColLonRadDeg ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], 
      mshCol, mshLon, mshRad );
    
    if ( (mshCol >= mod.colMin && mshCol <= mod.colMax) &&
          (mshLon >= mod.lonMin && mshLon <= mod.lonMax) &&
            (mshRad >= mod.radMin) ) {
                        
      kdres *set   = kd_nearest3 ( tree, msh.xmsh[i], msh.ymsh[i], 
        msh.zmsh[i] );    
      void  *ind_p = kd_res_item_data ( set );    
      int point    = * ( int * ) ind_p;
      
      double tap = taper ( mshCol, mshLon, mshRad, mod );
    
      msh.c11[i] = tap * mod.c11[point]    + msh.c11[i];
      msh.c12[i] = tap * mod.c12[point]    + msh.c12[i];
      msh.c13[i] = tap * mod.c13[point]    + msh.c13[i];
      msh.c22[i] = tap * mod.c22[point]    + msh.c22[i];
      msh.c23[i] = tap * mod.c23[point]    + msh.c23[i];
      msh.c33[i] = tap * mod.c33[point]    + msh.c33[i];
      msh.c44[i] = tap * mod.c44[point]    + msh.c44[i];
      msh.c55[i] = tap * mod.c55[point]    + msh.c55[i];
      msh.c66[i] = tap * mod.c66[point]    + msh.c66[i];
      msh.rho[i] = tap * mod.rhoMsh[point] + msh.rho[i];
      
    }    
  }
          
}

void Interpolator::exterpolator ( Mesh &msh, Exodus_file &exo, Model_file &mod )
{
  
  double c11, c12, c13, c14, c15, c16,c22, c23, c24, c25, c26, c33, c34, c35;
  double c36, c44, c45, c46, c55, c56, c66, rho;
  
  double testX;
  double testY;
  double testZ;
  
  Utilities util;
    
  // Create element connectivity map.  
  msh.getConnectivity ( exo.idexo );
  
  cout << "Extracting model.\n";  

  if ( mod.input_model_file_type == "SES3D" ) {
    
    // Create KDTree.
    cout << "Creating KDTree.\n";
    kdtree *tree = kd_create (3);
    int *dat     = new int [msh.num_nodes];
    for ( int i=0; i<msh.num_nodes; i++ ) {
      dat[i] = i;
      kd_insert3 ( tree, msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], &dat[i] );
    }                   

    cout << "Recovering values.\n";
    int l = 0;
    for ( int r=0; r<mod.col_deg.size(); r++ ) {
      for ( int i=0; i<mod.col_rad[r].size(); i++ ) {
        for ( int j=0; j<mod.lon_rad[r].size(); j++ ) {
          for ( int k=0; k<(mod.rad[r].size()-1); k++ ) {
                      
            // Test point (DEBUG).
            util.colLonRadRad2xyz ( mod.col_rad[r][i], mod.lon_rad[r][j], 
              mod.rad[r][k], testX, testY, testZ );
                          
            recover ( testX, testY, testZ, tree, msh,
              c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26,
              c33, c34, c35, c36, c44, c45, c46, c55, c56, c66,
              rho );
                      
            if ( mod.input_model_physics == "TTI" ) {
              mod.vsh[r][l] = sqrt (c44 / rho);
              mod.vsv[r][l] = sqrt (c55 / rho);
              mod.vpp[r][l] = sqrt (c22 / rho);
              mod.rho[r][l] = rho;
              l++;              
                            
            }
          }
        }
        cout << "CoLat loops left in this region: " << 
          ( mod.col_rad[r].size() - (i + 1) ) << "\xd" << std::flush;
      }  
    }
  }
   
  if ( mod.input_model_file_type == "SPECFEM" ) {
    // TODO Write specfem in.    
  }
  
}

double Interpolator::taper ( double &col, double &lon, double &rad, 
  Model_file &mod )
{
  
  Constants con;
  Utilities util;
  
  double tap;
  double dTaper = 500.;
  
  col = col * con.PI / con.o80;
  lon = lon * con.PI / con.o80;
    
  double lat       = util.col2Lat ( col, 'r' );  
  double modLatMin = util.col2Lat ( mod.colMin, 'd') * con.PI / con.o80;
  double modLatMax = util.col2Lat ( mod.colMax, 'd') * con.PI / con.o80;
  double modLonMin = mod.lonMin * con.PI / con.o80;
  double modLonMax = mod.lonMax * con.PI / con.o80;
      
  // Lat dist.
  double dLat = lat - modLatMin;
  double dsqt = pow ( sin(dLat / 2.), 2);
  double dphi = 2 * asin ( sqrt (dsqt) );
  
  double distLat = rad * dphi;
  
  // Lon dist.
  double dLon = lon - modLonMin;
         dsqt = cos (lat) * cos (lat) * pow ( sin(dLon / 2.), 2);
         dphi = 2 * asin ( sqrt (dsqt) );
         
  double distLon = rad * dphi;    
  
  // Max Lat.
  dLat = modLatMax - modLatMin;
  dsqt = pow ( sin(dLat / 2.), 2);
  dphi = 2 * asin ( sqrt (dsqt) );
  
  double distLatMax = rad * dphi - dTaper;
  
  // Max Lon.
  dLon = modLonMax - modLonMin;
  dsqt = cos (lat) * cos (lat) * pow (sin (dLon / 2.), 2);
  dphi = 2 * asin ( sqrt (dsqt) );
  
  double distLonMax = rad * dphi - dTaper;
  
  // Rad.
  double dRad = rad - mod.radMin;
  
  double tapLat = 1.;
  double tapLon = 1.;  
  double tapRad = 1.;
  if ( distLat <= dTaper ) {
    tapLat = distLat / dTaper;
  } else if ( distLat >= distLatMax ) {
    tapLat = (distLatMax + dTaper - distLat) / dTaper;
  }
  
  if ( distLon <= dTaper ) {
    tapLon = distLon / dTaper;    
  } else if ( distLon >= distLonMax ) {
    tapLon = (distLonMax + dTaper - distLon) / dTaper;
  }
  
  if ( dRad <= dTaper ) {
    tapRad = dRad / dTaper;
  }
      
  return tapLat * tapLon * tapRad;
    
}

void Interpolator::recover ( double &testX, double &testY, double &testZ, 
               kdtree *tree,
               Mesh &msh,
               double &c11, double &c12, double &c13, double &c14, double &c15, 
               double &c16, double &c22, double &c23, double &c24, double &c25, 
               double &c26, double &c33, double &c34, double &c35, double &c36, 
               double &c44, double &c45, double &c46, double &c55, double &c56, 
               double &c66, double &rho )
{
  
  Utilities util;
  
  /* orig* holds the originally requested points. The test* variables are used
  to recursively search in a random region around the selected point if we don't
  manage to find it on the first try */
  double origX = testX;
  double origY = testY;
  double origZ = testZ;

  // These variables hold the orignally requested col, lon, and lat.
  double col, lon, rad;
  double colPoint, lonPoint, radPoint;
  
  /* Here we keep track of whether we've found the variable, and if this is our
  first try. Assume we haven't found it at first.*/
  bool found=false;  
  bool first=true;
  
  while ( found == false ) {
        
    // Extract point from KDTree.
    kdres *set  = kd_nearest3 ( tree, testX, testY, testZ );
    void *ind_p = kd_res_item_data ( set );
    int point   = * ( int * ) ind_p;
    
    // Find the originally requested col, lon, and rad.
    if ( first == true ) {      
      util.xyz2ColLonRadRad ( origX, origY, origZ, col, lon, rad );        
      first = false;      
    }

    // Set up connectivity iterator.
    pair < multimap <int, vector <int> > :: iterator , multimap 
      <int, vector <int> > :: iterator > ext;

    // Extract iterator for current point.
    ext = msh.elemOrder.equal_range (point);

    // Loop over connecting elements (node indices contained in iterator).
    int nFound = 0;    
    double l1, l2, l3, l4;
    for ( multimap <int, vector <int> > :: iterator it=ext.first; 
      it!=ext.second; ++it ) {                    

      /* Convert to barycentric coordinates (l*). The vector second contains 
        all the indices of the nodes belonging to a element (4 for a tet). So
        we need to extract 12 values, 4 for each dimension */
      util.convertBary ( origX, origY, origZ,
        msh.xmsh[it->second[0]], msh.xmsh[it->second[1]], 
        msh.xmsh[it->second[2]], msh.xmsh[it->second[3]],
        msh.ymsh[it->second[0]], msh.ymsh[it->second[1]],
        msh.ymsh[it->second[2]], msh.ymsh[it->second[3]],
        msh.zmsh[it->second[0]], msh.zmsh[it->second[1]],
        msh.zmsh[it->second[2]], msh.zmsh[it->second[3]],
        l1, l2, l3, l4 ); 

      // If barycentric coordinates are all >= 0.
      if ( l1 >= 0 && l2 >= 0 && l3 >= 0 && l4 >= 0 ) {
      
        found = true;
      
        double c11p0 = msh.c11[it->second[0]];
        double c12p0 = msh.c12[it->second[0]];
        double c13p0 = msh.c13[it->second[0]];
        double c14p0 = msh.c14[it->second[0]];
        double c15p0 = msh.c15[it->second[0]];
        double c16p0 = msh.c16[it->second[0]];
        double c22p0 = msh.c22[it->second[0]];
        double c23p0 = msh.c23[it->second[0]];
        double c24p0 = msh.c24[it->second[0]];
        double c25p0 = msh.c25[it->second[0]];
        double c26p0 = msh.c26[it->second[0]];
        double c33p0 = msh.c33[it->second[0]];
        double c34p0 = msh.c34[it->second[0]];
        double c35p0 = msh.c35[it->second[0]];
        double c36p0 = msh.c36[it->second[0]];
        double c44p0 = msh.c44[it->second[0]];
        double c45p0 = msh.c45[it->second[0]];
        double c46p0 = msh.c46[it->second[0]];
        double c55p0 = msh.c55[it->second[0]];
        double c56p0 = msh.c56[it->second[0]];
        double c66p0 = msh.c66[it->second[0]]; 
        double rhop0 = msh.rho[it->second[0]];  
                            
        double c11p1 = msh.c11[it->second[1]];
        double c12p1 = msh.c12[it->second[1]];
        double c13p1 = msh.c13[it->second[1]];
        double c14p1 = msh.c14[it->second[1]];
        double c15p1 = msh.c15[it->second[1]];
        double c16p1 = msh.c16[it->second[1]];
        double c22p1 = msh.c22[it->second[1]];
        double c23p1 = msh.c23[it->second[1]];
        double c24p1 = msh.c24[it->second[1]];
        double c25p1 = msh.c25[it->second[1]];
        double c26p1 = msh.c26[it->second[1]];
        double c33p1 = msh.c33[it->second[1]];
        double c34p1 = msh.c34[it->second[1]];
        double c35p1 = msh.c35[it->second[1]];
        double c36p1 = msh.c36[it->second[1]];
        double c44p1 = msh.c44[it->second[1]];
        double c45p1 = msh.c45[it->second[1]];
        double c46p1 = msh.c46[it->second[1]];
        double c55p1 = msh.c55[it->second[1]];
        double c56p1 = msh.c56[it->second[1]];
        double c66p1 = msh.c66[it->second[1]];
        double rhop1 = msh.rho[it->second[1]];     
      
        double c11p2 = msh.c11[it->second[2]];
        double c12p2 = msh.c12[it->second[2]];
        double c13p2 = msh.c13[it->second[2]];
        double c14p2 = msh.c14[it->second[2]];
        double c15p2 = msh.c15[it->second[2]];
        double c16p2 = msh.c16[it->second[2]];
        double c22p2 = msh.c22[it->second[2]];
        double c23p2 = msh.c23[it->second[2]];
        double c24p2 = msh.c24[it->second[2]];
        double c25p2 = msh.c25[it->second[2]];
        double c26p2 = msh.c26[it->second[2]];
        double c33p2 = msh.c33[it->second[2]];
        double c34p2 = msh.c34[it->second[2]];
        double c35p2 = msh.c35[it->second[2]];
        double c36p2 = msh.c36[it->second[2]];
        double c44p2 = msh.c44[it->second[2]];
        double c45p2 = msh.c45[it->second[2]];
        double c46p2 = msh.c46[it->second[2]];
        double c55p2 = msh.c55[it->second[2]];
        double c56p2 = msh.c56[it->second[2]];
        double c66p2 = msh.c66[it->second[2]];    
        double rhop2 = msh.rho[it->second[2]];     
      
        double c11p3 = msh.c11[it->second[3]];
        double c12p3 = msh.c12[it->second[3]];
        double c13p3 = msh.c13[it->second[3]];
        double c14p3 = msh.c14[it->second[3]];
        double c15p3 = msh.c15[it->second[3]];
        double c16p3 = msh.c16[it->second[3]];
        double c22p3 = msh.c22[it->second[3]];
        double c23p3 = msh.c23[it->second[3]];
        double c24p3 = msh.c24[it->second[3]];
        double c25p3 = msh.c25[it->second[3]];
        double c26p3 = msh.c26[it->second[3]];
        double c33p3 = msh.c33[it->second[3]];
        double c34p3 = msh.c34[it->second[3]];
        double c35p3 = msh.c35[it->second[3]];
        double c36p3 = msh.c36[it->second[3]];
        double c44p3 = msh.c44[it->second[3]];
        double c45p3 = msh.c45[it->second[3]];
        double c46p3 = msh.c46[it->second[3]];
        double c55p3 = msh.c55[it->second[3]];
        double c56p3 = msh.c56[it->second[3]];
        double c66p3 = msh.c66[it->second[3]];    
        double rhop3 = msh.rho[it->second[3]];     
      
        c11 = l1 * c11p0 + l2 * c11p1 + l3 * c11p2 + l4 * c11p3;    
        c12 = l1 * c12p0 + l2 * c12p1 + l3 * c12p2 + l4 * c12p3;    
        c13 = l1 * c13p0 + l2 * c13p1 + l3 * c13p2 + l4 * c13p3;    
        c14 = l1 * c14p0 + l2 * c14p1 + l3 * c14p2 + l4 * c14p3;    
        c15 = l1 * c15p0 + l2 * c15p1 + l3 * c15p2 + l4 * c15p3;    
        c16 = l1 * c16p0 + l2 * c16p1 + l3 * c16p2 + l4 * c16p3;    
        c22 = l1 * c22p0 + l2 * c22p1 + l3 * c22p2 + l4 * c22p3;    
        c23 = l1 * c23p0 + l2 * c23p1 + l3 * c23p2 + l4 * c23p3;    
        c24 = l1 * c24p0 + l2 * c24p1 + l3 * c24p2 + l4 * c24p3;    
        c25 = l1 * c25p0 + l2 * c25p1 + l3 * c25p2 + l4 * c25p3;    
        c26 = l1 * c26p0 + l2 * c26p1 + l3 * c26p2 + l4 * c26p3;    
        c33 = l1 * c33p0 + l2 * c33p1 + l3 * c33p2 + l4 * c33p3;    
        c34 = l1 * c34p0 + l2 * c34p1 + l3 * c34p2 + l4 * c34p3;    
        c35 = l1 * c35p0 + l2 * c35p1 + l3 * c35p2 + l4 * c35p3;    
        c36 = l1 * c36p0 + l2 * c36p1 + l3 * c36p2 + l4 * c36p3;    
        c44 = l1 * c44p0 + l2 * c44p1 + l3 * c44p2 + l4 * c44p3;    
        c45 = l1 * c45p0 + l2 * c45p1 + l3 * c45p2 + l4 * c45p3;    
        c46 = l1 * c46p0 + l2 * c46p1 + l3 * c46p2 + l4 * c46p3;    
        c55 = l1 * c55p0 + l2 * c55p1 + l3 * c55p2 + l4 * c55p3;    
        c56 = l1 * c56p0 + l2 * c56p1 + l3 * c56p2 + l4 * c56p3;    
        c66 = l1 * c66p0 + l2 * c66p1 + l3 * c66p2 + l4 * c66p3;  
        rho = l1 * rhop0 + l2 * rhop1 + l3 * rhop2 + l4 * rhop3; 
                
        break;              
       }                     
    }
    
    /* If we haven't found the point on the first try, need to do some magical
    stuff */
    if ( found == false ) {                
      
      /* For col and lon, randomly choose which direction to look. This might
      be able to be switched to a more direction search (i.e. we look in the
      direction which is towards the requested point), but this is certainly
      a bit more robust. For radius, we do a directional search (i.e. search 
      only down if we are higher in radius). This is because radial 
      discretization is much finer than lat/lon. */
      int signC = rand () % 2 ? 1 : -1;  
      int signL = rand () % 2 ? 1 : -1;
      int signR = radPoint < rad ? 1 : -1;
      
      // TODO Make this average edge length a variable.
      /* The theta discretization is controlled by the average edge length of an
      element in the lat/lon direction */      
      double dTheta = 85 / rad;
      
      /* Allow the search radius to range from 0 to 1 times some values */
      double randC = (rand () % 100) / 100.;
      double randL = (rand () % 100) / 100.;      
      double randR = (rand () % 100) / 100.;
      
      /* Col and Lon are allowed to vary between their original values (0) and
      one edge length away. This seems to work. Radius is allowed to vary by
      1 km. This also seems to work, although I think it's a bit sketchier. 
      Could add a parameter to the radius search to make this variable, 
      dependent on depth. It does work quite well now though */
      double colTest = col + ( signC * randC * dTheta );
      double lonTest = lon + ( signL * randL * dTheta );
      double radTest = rad + ( signR * randR );
      
      /* Create a new testX, Y, and Z, point for the recursive search. */
      util.colLonRadRad2xyz ( colTest, lonTest, radTest, testX, testY, testZ );
      
    }    
  }  
}  