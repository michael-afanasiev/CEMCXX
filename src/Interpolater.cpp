#include <iostream>
#include <cmath>
#include <assert.h>
#include <ctype.h>

#include "kdtree.h"
#include "classes.hpp"

using namespace std;

void Interpolator::interpolate ( Mesh &msh, Model_file &mod ) 
{
  
  kdtree *tree = kd_create (3);
  Utilities util;
  
  int *dat = new int [mod.num_p];
  for ( int i=0; i<mod.num_p; i++ ) {
    dat[i] = i;
    kd_insert3 ( tree, mod.x[i], mod.y[i], mod.z[i], &dat[i] );
    cout << "Creating tree. " << mod.num_p - (i + 1) << " nodes left.\xd";
  }      
  
  cout << "\n";
  
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
    
      msh.c11[i] = mod.c11[point] + msh.c11[i];
      msh.c12[i] = mod.c12[point] + msh.c12[i];
      msh.c13[i] = mod.c13[point] + msh.c13[i];
      msh.c22[i] = mod.c22[point] + msh.c22[i];
      msh.c23[i] = mod.c23[point] + msh.c23[i];
      msh.c33[i] = mod.c33[point] + msh.c33[i];
      msh.c44[i] = mod.c44[point] + msh.c44[i];
      msh.c55[i] = mod.c55[point] + msh.c55[i];
      msh.c66[i] = mod.c66[point] + msh.c66[i];
      
    }
    
    cout << "Interpolating. " << msh.num_nodes - (i + 1) << " nodes left.\xd";
    
  }
        
}

void Interpolator::exterpolator ( Mesh &msh, Exodus_file &exo, Model_file &mod )
{
  
  kdtree *tree = kd_create (3);
  Utilities util;
  
  int *dat = new int [mod.num_p];
  for ( int i=0; i<mod.num_p; i++ ) {
    dat[i] = i;
    kd_insert3 ( tree, mod.x[i], mod.y[i], mod.z[i], &dat[i] );
    cout << "Creating tree. " << mod.num_p - (i + 1) << " nodes left.\xd";
  }      
  
  cout << "\n";
  
  // Create element connectivity map.  
  msh.getConnectivity ( exo.idexo );
  
  double testX = 1000.;
  double testY = 1000.;
  double testZ = 1000.;
  
  kdres *set  = kd_nearest3 ( tree, testX, testY, testZ );
  void *ind_p = kd_res_item_data ( set );
  int point   = * ( int * ) ind_p;
  
  cout << "POINT IS: " << point;
  
  pair < multimap <int, vector <int> > :: iterator , multimap 
    <int, vector <int> > :: iterator > ext;
  
  ext = msh.elemOrder.equal_range (point);
  
  for ( multimap <int, vector <int> > :: iterator it=ext.first; 
    it!=ext.second; ++it ) {                    
    cout << it->second[0] << "\n" << it->second[1] << "\n" << it-> second[2] << "\n"
      << it -> second[3] << "\n";
    cin.get ();
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