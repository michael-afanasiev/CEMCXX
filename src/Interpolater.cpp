#include <iostream>
#include <cmath>
#include <assert.h>
#include <ctype.h>
#include <time.h>

#include "kdtree.h"
#include "classes.hpp"

using namespace std;

void Interpolator::interpolateCrust ( Mesh &msh, Discontinuity &dis )
{
  
  Utilities utl;
  
  cout << "Adding crust.\n";
  
  if ( msh.radMin > 6271 )
  {
    for ( int i=0; i<msh.num_nodes; i++ )
    {
    
      double col, lon, rad;
    
      bool inCrust = false;
    
      utl.xyz2ColLonRadDeg ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], 
        col, lon, rad);
      
      dis.lookCrust ( msh, col, lon, rad, i, inCrust );    
    }
  }
}

void Interpolator::interpolateTopo ( Mesh &msh, Discontinuity &dis )
{
  
  Utilities utl;
  Constants con;
  
  cout << "Adding topography." << endl;
  
  if ( msh.radMin > (6271.) )
  {
    for ( int i=0; i<msh.num_nodes; i++ )
    {
      double col, lon, rad;
        
      utl.xyz2ColLonRadDeg ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], 
        col, lon, rad);
        
      dis.lookTopo ( msh, col, lon, rad, i );
    }
  }
}

void Interpolator::findNodes ( Mesh &msh, Model_file &mod, ofstream &myfile )
{
  
  Utilities utl;
  Constants con;
  
  cout << "Finding nodes to refine." << endl;
  
  int outCount = 0; 
  int inCount  = 0;
  
  std::vector <int> node;
  for ( std::vector < vector <int> > :: iterator out=msh.refineElemConn.begin();
    out!=msh.refineElemConn.end(); ++out )       
  {
    
    for ( std::vector <int> :: iterator in=out->begin(); in!=out->end(); ++in )
    {
      
      // Define local variables.
      double mshColRot, mshLonRot, mshRadRot; // Sph. Coord. in rotated domain.
      double xRot,      yRot,      zRot;      // Cart. Coord in rotated domain.
  
      /* Rotate from physical domain ( in exodus file ) to simulation domain
      ( in SES3D file ). Grab points in simulation domain */
      utl.rotateBackward ( msh.xmsh[*in-1], msh.ymsh[*in-1], msh.zmsh[*in-1], 
        xRot, yRot, zRot, mod );
        
      utl.xyz2ColLonRadDeg ( xRot, yRot, zRot, mshColRot, mshLonRot, 
        mshRadRot );
    
      // Handle special cases of longitude axis wrapping.
      if ( mod.wrapAround == true && mshLonRot < 0. )
      {
        mshLonRot += 360;   
      }
      
      outCount++;
      if ( (mshColRot >= mod.colMin && mshColRot <= mod.colMax) &&
           (mshLonRot >= mod.lonMin && mshLonRot <= mod.lonMax) &&
           (mshRadRot >= mod.radMin) ) 
      {       
        inCount++;
      
        if ( inCount == msh.num_node_per_elem )
        {
          node.push_back ( *(in-3) );
          node.push_back ( *(in-2) );
          node.push_back ( *(in-1) );
          node.push_back ( *(in-0) );  
        }

      }          
      
      if ( (outCount % msh.num_node_per_elem) == 0 )
      {
        inCount = 0;
      }
            
    }                        
    
  }
  

  if ( node.size() != 0 )
  {
    msh.refineElemConn.push_back (node);
    cout << "Pushed" << endl;
  }
  
}

void Interpolator::interpolate ( Mesh &msh, Model_file &mod, Discontinuity 
  &dis ) 
{

  Attenuation atn;
  Utilities utl;
  Constants con;
  Mod1d     bm;
    
  cout << "Interpolating.\n";
    
  // Loop over every node point in the exodus file.
#pragma omp parallel for
  for ( int i=0; i<msh.num_nodes; i++ ) {        

    // Define local variables.
    double mshColRot, mshLonRot, mshRadRot; // Sph. Coord. in rotated domain.
    double mshColPys, mshLonPys, mshRadPys; // Sph. Coord. in phsyical domain.
    double xRot,      yRot,      zRot;      // Cart. Coord in rotated domain.
    
    // At each node point, assume we are not in the crust.
    bool inCrust = false;
    
    /* Rotate from physical domain ( in exodus file ) to simulation domain
    ( in SES3D file ). Grab points in simulation domain */
    utl.rotateBackward ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], xRot, 
      yRot, zRot, mod );
    
    /* Get sph. coord. for both physical domain ( for discontinutiy search )
      and rotated domain ( to ensure we are only interpolating points onto the 
      mesh which are coincident to the interpolated file ). */
    utl.xyz2ColLonRadDeg ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], 
      mshColPys, mshLonPys, mshRadPys );
    utl.xyz2ColLonRadDeg ( xRot, yRot, zRot, mshColRot, mshLonRot, mshRadRot );
    
    // TODO this is fucked up
    // Handle special cases of longitude axis wrapping.
    // if ( mod.wrapAround == true && mshLonPys < 0. )
    //   mshLonPys += 360;        
    // if ( mod.wrapAround == true && mshLonRot < 0. )
    //   mshLonRot += 360;    
      
    // Check for discontinuity conditinos.
    utl.checkRegion ( msh, mshRadPys );
    dis.lookCrust   ( msh, mshColPys, mshLonPys, mshRadPys, i, inCrust );

    /* If the rotated coordinates are within the simulation domain, go ahead and
    interpolate. */
    if ( (mshColRot >= mod.colMin && mshColRot <= mod.colMax) &&
         (mshLonRot >= mod.lonMin && mshLonRot <= mod.lonMax) &&
         (mshRadRot >= mod.radMin) ) 
    {                                             
                  
      /* Find the xyz values which are closest to the rotated point. Index of 
      the xyz values are stored in 'point'. */       
      kdres *set   = kd_nearest3 ( mod.tree, xRot, yRot, zRot );    
      void  *ind_p = kd_res_item_data ( set );    
      int point    = * ( int * ) ind_p;
                               
      // Check taper condition based on distance from edge of rotated model.
      double tap = taper ( mshColRot, mshLonRot, mshRadRot, mod );
      
      // Adjust radius and ensure we don't grab from across boundaries.
      utl.checkRegion ( msh, mshRadRot );
      
      // Get 1d background values.
      double vs1d, vp1d, rho1d;
      bm.eumod                   ( mshRadRot, vs1d, vp1d, rho1d );     
      double qvCor = atn.correct ( atn.qModelName, mshRadRot );
      
      // FIXME 
      // qvCor == 1.;
                   
      // TTI.
      double vshExo = sqrt ( msh.c44[i] / msh.rho[i] );
      double vsvExo = sqrt ( msh.c55[i] / msh.rho[i] );
      double vppExo = sqrt ( msh.c22[i] / msh.rho[i] );
      double rhoExo = msh.rho[i];
                  
      double vshMod = mod.vshUnwrap[point];
      double vsvMod = mod.vsvUnwrap[point];
      double vppMod = mod.vppUnwrap[point];
      double rhoMod = mod.rhoUnwrap[point]; 

      double rhoModCor;
      double vshModCor;
      double vsvModCor;
      double vppModCor;

      if ( mod.kernel == false )
      {        
        rhoModCor = rhoMod;
        vshModCor = vshMod * qvCor;
        vsvModCor = vsvMod * qvCor;
        vppModCor = vppMod * qvCor;
      }
      else
      {
        rhoModCor = ( rhoMod + rho1d );      
        vshModCor = ( vshMod + vs1d ) * qvCor;
        vsvModCor = ( vsvMod + vs1d ) * qvCor;
        vppModCor = ( vppMod + vp1d ) * qvCor;
      }
            
      double rhoUse = ( 1 - tap ) * rhoExo + tap * rhoModCor;
      double vshUse = ( 1 - tap ) * vshExo + tap * vshModCor;
      double vsvUse = ( 1 - tap ) * vsvExo + tap * vsvModCor;
      double vppUse = ( 1 - tap ) * vppExo + tap * vppModCor;
            
      double N = rhoUse * vshUse * vshUse;
      double L = rhoUse * vsvUse * vsvUse;
      double A = rhoUse * vppUse * vppUse;
      
      double C = A;
      double F = A - 2 * L;
      double S = A - 2 * N;

      if ( (inCrust == false) || 
           ((mod.overwriteCrust == true) && (dis.inCrust == true)) ) 
      {
        cout << "hi" << endl;
        msh.c11[i] = C;
        msh.c22[i] = A;
        msh.c33[i] = A;
        msh.c12[i] = F;
        msh.c13[i] = F;
        msh.c23[i] = S;
        msh.c44[i] = N;
        msh.c55[i] = L;
        msh.c66[i] = L;
        msh.rho[i] = rhoUse;              
        
      }                  
    }          
  }          
}

double Interpolator::taper ( double &col, double &lon, double &rad, 
  Model_file &mod )
{
  
  Constants con;
  Utilities util;
  
  double dTaper    = 2000.;
  double dTaperRad = 50.;
  
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
  
  if ( dRad <= dTaperRad ) {
    tapRad = dRad / dTaperRad;
  }
      
  return tapLat * tapLon * tapRad;
    
}

int Interpolator::recover ( double &testX, double &testY, double &testZ, 
               Mesh &msh, double &c11, double &c12, double &c13, double &c14, 
               double &c15, double &c16, double &c22, double &c23, double &c24, 
               double &c25, double &c26, double &c33, double &c34, double &c35, 
               double &c36, double &c44, double &c45, double &c46, double &c55, 
               double &c56, double &c66, double &rho, char mode, 
               bool &internalFound )
{
  
  Utilities util;
  Constants con;
  
  internalFound = false;
  
  kdtree *tree = msh.tree;
  kdres  *set;
  void   *ind_p;
  
  /* orig* holds the originally requested points. The test* variables are used
  to recursively search in a random region around the selected point if we don't
  manage to find it on the first try */
  double origX = testX;
  double origY = testY;
  double origZ = testZ;

  // These variables hold the orignally requested col, lon, and lat.
  double col, lon, rad;
  double colPoint, lonPoint, radPoint;
  double colClose, lonClose, radClose;
  
  /* Here we keep track of whether we've found the variable, and if this is our
  first try. Assume we haven't found it at first.*/
  // std::vector < int > repeater;
  // repeater.reserve ( 100 );
  bool found   = false;  
  int  count   = 0;
  int  nodeNum = 0;
  int  point;
  
#ifdef VISUAL_DEBUG
  ofstream myfile;
  myfile.open ( "Bary.txt", ios::out );
#endif
                                
  // Find node closest to point.
  set   = kd_nearest3 ( tree, testX, testY, testZ );
  ind_p = kd_res_item_data ( set );
  point = * ( int * ) ind_p;
  
  kd_res_free ( set );
            
  // Find the originally requested col, lon, and rad.
  util.xyz2ColLonRadRad ( origX, origY, origZ, col, lon, rad );        
  
  // Find ColLonRad of original node.
  util.xyz2ColLonRadRad ( msh.xmsh[point], msh.ymsh[point], msh.zmsh[point], 
    colClose, lonClose, radClose );

  // Repeat until enclosing tet is found.
  while ( found == false ) {    
  
    // Set up connectivity iterator.
    pair < multimap <int, vector <int> > :: iterator , multimap 
      <int, vector <int> > :: iterator > ext;

    // Extract iterator for current point.
    ext = msh.elemOrder.equal_range (point);
    
    // Extract node number if looking for it.
    if ( mode == 'e' )
    {
      nodeNum = msh.node_num_map[point];
    }
    
    // Loop over connecting elements (node indices contained in iterator).
    double l1, l2, l3, l4;
    for ( multimap <int, vector <int> > :: iterator it=ext.first; 
      it!=ext.second; ++it ) 
      {                    

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
          
#ifdef VISUAL_DEBUG
        double col1, lon1, rad1;
        double col2, lon2, rad2;
        double col3, lon3, rad3;
        double col4, lon4, rad4;
              
        cout << l1 << " " << l2 << " " << l3 << " " << l4 << endl;
        
        util.xyz2ColLonRadDeg ( msh.xmsh[it->second[0]],  
                                msh.ymsh[it->second[0]],
                                msh.zmsh[it->second[0]],
                                col1, lon1, rad1 );
        util.xyz2ColLonRadDeg ( msh.xmsh[it->second[1]],  
                                msh.ymsh[it->second[1]],
                                msh.zmsh[it->second[1]],
                                col2, lon2, rad2 );
        util.xyz2ColLonRadDeg ( msh.xmsh[it->second[2]],  
                                msh.ymsh[it->second[2]],
                                msh.zmsh[it->second[2]],
                                col3, lon3, rad3 );
        util.xyz2ColLonRadDeg ( msh.xmsh[it->second[3]],  
                                msh.ymsh[it->second[3]],
                                msh.zmsh[it->second[3]],
                                col4, lon4, rad4 );
                               
                               
        cout << col1 << " " << lon1 << " " << rad1 << endl;                      
        cout << col2 << " " << lon2 << " " << rad2 << endl;                      
        cout << col3 << " " << lon3 << " " << rad3 << endl;                      
        cout << col4 << " " << lon4 << " " << rad4 << endl;                      
        cout << col * con.o80 / con.PI  << " " << lon * con.o80 / con.PI << 
          " " << rad  << endl;
                                              
        myfile << 0 << " " << 0 << " " << origX << " " << origY << " " 
          << origZ << endl;
        myfile << 0 << " " << 0 << " " << msh.xmsh[it->second[0]] << " " << 
          msh.ymsh[it->second[0]] << " " << msh.zmsh[it->second[0]] << endl;
        
        myfile << 0 << " " << 0 << " " << msh.xmsh[it->second[1]] << " " << 
          msh.ymsh[it->second[1]] << " " << msh.zmsh[it->second[1]] << endl;
        
        myfile << 0 << " " << 0 << " " << msh.xmsh[it->second[2]] << " " << 
          msh.ymsh[it->second[2]] << " " << msh.zmsh[it->second[2]] << endl;
        
        myfile << 0 << " " << 0 << " " << msh.xmsh[it->second[3]] << " " << 
          msh.ymsh[it->second[3]] << " " << msh.zmsh[it->second[3]] << endl;
#endif
        
      // If barycentric coordinates are all >= 0.
      if ( l1 >= 0 && l2 >= 0 && l3 >= 0 && l4 >= 0 ) {
      
        found         = true;
        internalFound = true;
        
        if ( mode == 'p' ) 
        {
      
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
    }
    
    /* If we haven't found the point on the first try, need to do some magical
    stuff */
    if ( found == false ) 
    {   
      
      util.xyz2ColLonRadRad ( testX, testY, testZ, colPoint, lonPoint, 
        radPoint );             
      
      count++;
      // cout << count << endl;
      
      /* For col and lon, randomly choose which direction to look. This might
      be able to be switched to a more direction search (i.e. we look in the
      direction which is towards the requested point), but this is certainly
      a bit more robust. For radius, we do a directional search (i.e. search 
      only down if we are higher in radius). This is because radial 
      discretization is much finer than lat/lon. */
      int signC = rand () % 2 ? 1 : -1;  
      int signL = rand () % 2 ? 1 : -1;
      int signR = radClose < rad ? 1 : -1;
      
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
      if ( count < 50 ) 
      {
        double colTest = col + ( signC * randC * dTheta );
        double lonTest = lon + ( signL * randL * dTheta );
        double radTest = rad + ( signR * randR );
        
        /* Create a new testX, Y, and Z, point for the recursive search. */
        util.colLonRadRad2xyz ( colTest, lonTest, radTest, testX, testY, 
          testZ );          
          
        // Extract point from KDTree.
        set  = kd_nearest3 ( tree, testX, testY, testZ );
        ind_p = kd_res_item_data ( set );
        point       = * ( int * ) ind_p;  
        
        kd_res_free ( set );
        
        
      }
      else if ( count < 100 )
      {
        testX = origX + ( signC * randC * 85 );
        testY = origY + ( signL * randL * 85 );
        testZ = origZ + ( signR * randR * 50 );
        
        // Extract point from KDTree.
        set  = kd_nearest3 ( tree, testX, testY, testZ );
        ind_p = kd_res_item_data ( set );
        point       = * ( int * ) ind_p;  
        
        kd_res_free ( set );
                
      }
      
      if ( mode == 's' || count >= 100 )
      {
        // cout << "Not in this file bru. Taking closest point." << endl;
        // cout << "Col, lon, rad: " << col * con.o80 / con.PI << 
        //   " " << lon * con.o80 / con.PI << " " << rad << endl;
        // cout << origX << " " << origY << " " << origZ << endl;
        
        // Extract point from KDTree.
        set  = kd_nearest3 ( tree, origX, origY, origZ );
        ind_p = kd_res_item_data ( set );
        point       = * ( int * ) ind_p;  
        
        kd_res_free ( set );        
        
        c11 = msh.c11[point];
        c12 = msh.c12[point];
        c13 = msh.c13[point];
        c14 = msh.c14[point];
        c15 = msh.c15[point];
        c16 = msh.c16[point];
        c22 = msh.c22[point];
        c23 = msh.c23[point];
        c24 = msh.c24[point];
        c25 = msh.c25[point];
        c26 = msh.c26[point];
        c33 = msh.c33[point];
        c34 = msh.c34[point];
        c35 = msh.c35[point];
        c36 = msh.c36[point];
        c44 = msh.c44[point];
        c45 = msh.c45[point];
        c46 = msh.c46[point];
        c55 = msh.c55[point];
        c56 = msh.c56[point];
        c66 = msh.c66[point];    
        rho = msh.rho[point];    
        
        if ( rho < 1 )
        {
          cout << rho << " " << rad << endl;
        }
        
        
        found = true;      
        break;         
      }
                          
      if ( count >= 100 ) 
      {
        cout << "Looping forever. And ever. " << 
          "Boring! I'm outta here." << endl;
        
        cout << "I suggest you look into 'Interpolater.cpp' and adjust the " <<
          "search limits." << endl;
        
        cout.precision(15);
        cout << "Col, lon, rad: " << col * con.o80 / con.PI << 
          " " << lon * con.o80 / con.PI << " " << rad << endl;
        cout << "Col, lon, rad: " << col << 
          " " << lon  << " " << rad << endl;
        
#ifdef VISUAL_DEBUG
        myfile.close();
#endif
        
        exit ( EXIT_FAILURE );
      }
    }    
  }  
  
#ifdef VISUAL_DEBUG
  myfile.close();
#endif
          
  if ( mode == 'e' ) 
  {
    return nodeNum;
  } 
  else 
  { 
    // if ( count > 0 )
    //   {
    //   cout << count << endl;
    // }
    return 0;
  }
  
  kd_res_free ( set );
  
}  