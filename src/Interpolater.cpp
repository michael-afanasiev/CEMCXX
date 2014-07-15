#include <iostream>

#include <cmath>
#include <assert.h>
#include <ctype.h>
#include <time.h>

#include "kdtree.h"
#include "classes.hpp"

using namespace std;

void Interpolator::interpolateCrust ( Mesh &msh, Discontinuity &dis, 
  Model_file &mod )
{
  
  Utilities utl;
  
  std::cout << "Adding crust." << std::flush << std::endl;
  
#pragma omp parallel for
  for ( int i=0; i<msh.num_nodes; i++ )
  {
    double col, lon, rad, upTap, downTap;
    bool smoothCrust;
  
    bool inCrust = false;
  
    utl.xyz2ColLonRadDeg ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], 
      col, lon, rad);
    
    dis.lookCrust ( msh, col, lon, rad, i, inCrust, smoothCrust, upTap, 
      downTap, mod );     
  }
}

void Interpolator::interpolateTopo ( Mesh &msh, Discontinuity &dis )
{
  
  Utilities utl;
  
  std::cout << "Adding topography. " << std::flush << std::endl;

#pragma omp parallel for
  for ( int i=0; i<msh.num_nodes; i++ )
  {
    double col, lon, rad;
      
    utl.xyz2ColLonRadDeg ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], 
      col, lon, rad);
    dis.lookTopo ( msh, col, lon, rad, i );
  }
}

void Interpolator::interpolate ( Mesh &msh, Model_file &mod, Discontinuity 
  &dis ) 
{

  Attenuation atn;
  Utilities utl;
  Constants con;
  Mod1d     bm;
    
  std::cout << "Interpolating." << std::flush << std::endl; 

  // Loop over every node point in the exodus file.
#pragma omp parallel for
  for ( int i=0; i<msh.num_nodes; i++ ) {        

    kdres *set1;
    kdres *set2;
    void  *ind_p;
    int   point;

    // Define local variables.
    double mshColRot, mshLonRot, mshRadRot; // Sph. Coord. in rotated domain.
    double mshColPys, mshLonPys, mshRadPys; // Sph. Coord. in phsyical domain.
    double xRot,      yRot,      zRot;      // Cart. Coord in rotated domain.
    
    // At each node point, assume we are not in the crust.
    bool inCrust = false;
    
    /* Rotate from physical domain ( in exodus file ) to simulation domain
     * ( in SES3D file ). Grab points in simulation domain. Does not rotate if
     * we don't need to. */
    utl.rotateBackward ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], xRot, 
      yRot, zRot, mod );
    
    /* Get sph. coord. for both physical domain ( for discontinutiy search )
     * and rotated domain ( to ensure we are only interpolating points onto the 
     * mesh which are coincident to the interpolated file ). */
    utl.xyz2ColLonRadDeg ( msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], 
      mshColPys, mshLonPys, mshRadPys );
    utl.xyz2ColLonRadDeg ( xRot, yRot, zRot, mshColRot, mshLonRot, mshRadRot );
  
    // This should solve the wrap around issue I think.
    if ( mod.lonMax > 180. && mshLonPys <= 0. )
      mshLonPys += 360;
    if ( mod.lonMax > 180. && mshLonRot <= 0. )
      mshLonRot += 360;    
    if ( msh.lonMin < (-1 * con.PI / 2) && msh.lonMax > (con.PI / 2) )
      msh.lonMax = -1 * con.PI / 2;

    // Check for discontinuity conditions.
    double upTap, downTap;
    bool smoothCrust = false;
    utl.checkRegion ( msh, mshRadPys );
    dis.lookCrust   ( msh, mshColPys, mshLonPys, mshRadPys, i, inCrust, 
                      smoothCrust, upTap, downTap, mod );

    /* If the rotated coordinates are within the simulation domain, go ahead and
    interpolate. */
    if ( (mshColRot >= mod.colMin && mshColRot <= mod.colMax) &&
         (mshLonRot >= mod.lonMin && mshLonRot <= mod.lonMax) &&
         (mshRadRot >= mod.radMin) ) 
    {                                             
                  
      /* Find the xyz values which are closest to the rotated point. Index of 
      the xyz values are stored in 'point'. */       

      // Adjust radius and ensure we don't grab from across boundaries.
      utl.checkRegion ( msh, mshRadRot );
      
      /* This section is necessary due to the possibility of crossing region 
       * boundaries (specifically in a ses3d file). If we do cross boundaries
       * with a mesh file, two trees are automatically generated. Here we 
       * decide which one to extract from. */
      
      /* This parameter selects one of the two kdtrees, based on the actual
       * values of the mesh radius. It's set to -1 first as a failsafe. */
      int kdRegionExtract = -1;
      for ( int r=0; r<mod.kdRegions.size(); r++ )
      {

        /* Check which model region we are in by looping through them all
         * and seeing which radial chunk we're in. KDregions maps the actual
         * region number to the right kdTree. */
        // IF WE GET DIAMONDING, ADD THE -1 & +1 BACK HERE. 
        if ( mshRadRot >= ( mod.minRadReg[mod.kdRegions[r]] ) && 
             mshRadRot <= ( mod.maxRadReg[mod.kdRegions[r]] ) )
        {
          kdRegionExtract = r;
        }

      }

      /* If we're in 'no' region, that's obviously a problem. Exit. */
      if ( kdRegionExtract == -1 )
      {
        std::cout << "__FATAL__" << std::flush << std::endl;
        exit ( EXIT_FAILURE );
      }

      /* Extract from the proper kdtree, depending on the value extracted from
       * kdRegions above. */
      if ( kdRegionExtract == 0 )
      {
        set1  = kd_nearest3 ( mod.tree1, xRot, yRot, zRot );   
        ind_p = kd_res_item_data ( set1 );   
        point = * ( int * ) ind_p; 
        kd_res_free ( set1 );
      }

      if ( kdRegionExtract == 1 )
      {
        set2  = kd_nearest3 ( mod.tree2, xRot, yRot, zRot );   
        ind_p = kd_res_item_data ( set2 );   
        point = * ( int * ) ind_p; 
        kd_res_free ( set2 );
      }

      // Check taper condition based on distance from edge of rotated model.
      double tap = taper ( mshColRot, mshLonRot, mshRadRot, mod );

      // Get 1d background values.
      double vs1d, vp1d, rho1d;
      bm.eumod                   ( mshRadRot, vs1d, vp1d, rho1d );     
      
      // Fix attenuation.
      double qvCor = atn.correct ( atn.qModelName, mshRadRot );
      
      /* Project to different symmetry systems here */
      double N      = 0; 
      double L      = 0; 
      double A      = 0; 
      double C      = 0;
      double F      = 0; 
      double S      = 0; 
      double rhoUse = 0;                   
      if ( mod.input_model_physics == "TTI" )
      {  
        // Build according to symmetry system.
        double vshExo = sqrt ( msh.c44[i] / msh.rho[i] );
        double vsvExo = sqrt ( msh.c55[i] / msh.rho[i] );
        double vppExo = sqrt ( msh.c22[i] / msh.rho[i] );
        double rhoExo = msh.rho[i];
                    
        // These guys hold the model perturbations at a point.
        double vshMod = mod.vshUnwrap[point];
        double vsvMod = mod.vsvUnwrap[point];
        double vppMod = mod.vppUnwrap[point];
        double rhoMod = mod.rhoUnwrap[point]; 
        
        // These parameters hold values after correction to 1s for attenuation.
        double rhoModCor;
        double vshModCor;
        double vsvModCor;
        double vppModCor;
        
        /* Here we either a) replace the model, b) add the model as a 
         * perturbation to a 1d background, or c) add the model to the 3d mesh.
         * There's a small complication: attenuation. For the first two cases, 
         * it's fine. In the 3rd case, I add the model to a 1d background, 
         * compute the corrected value at 1s, and then subtract the 1d 
         * background. */
        /* TODO Check to see if we can add attenuation directly to the 3d model.
         * for instance, could probably multiple vshExo by 1/qvCor, add the
         * perturbation, and then multiply again */        
        if ( mod.kernel1d == false && mod.kernel3d == false )
        {        
          rhoModCor = rhoMod;
          vshModCor = vshMod * qvCor;
          vsvModCor = vsvMod * qvCor;
          vppModCor = vppMod * qvCor;
        }
        else if ( mod.kernel1d == true )
        {
          rhoModCor = ( rhoMod + rho1d );      
          vshModCor = ( vshMod + vs1d ) * qvCor;
          vsvModCor = ( vsvMod + vs1d ) * qvCor;
          vppModCor = ( vppMod + vp1d ) * qvCor;
        }
        else if ( mod.kernel3d == true )
        {
          double rhoModAtn = ( rhoMod + rho1d );
          double vshModAtn = ( vshMod + vs1d ) * qvCor;
          double vsvModAtn = ( vsvMod + vs1d ) * qvCor;
          double vppModAtn = ( vppMod + vp1d ) * qvCor;

          rhoModCor = rhoModAtn - rho1d + rhoExo;
          vshModCor = vshModAtn - vs1d  + vshExo;
          vsvModCor = vsvModAtn - vs1d  + vsvExo;
          vppModCor = vppModAtn - vp1d  + vppExo;
        }
              
        // Edge tapering.
        double vshUse = ( 1 - tap ) * vshExo + tap * vshModCor;
        double vsvUse = ( 1 - tap ) * vsvExo + tap * vsvModCor;
        double vppUse = ( 1 - tap ) * vppExo + tap * vppModCor;
        rhoUse        = ( 1 - tap ) * rhoExo + tap * rhoModCor;
         
        // TTI parameters to add.
        N = rhoUse * vshUse * vshUse;
        L = rhoUse * vsvUse * vsvUse;
        A = rhoUse * vppUse * vppUse;        
        F = A - 2 * L;
        S = A - 2 * N;
        C = A;
        
      }

      // In this case, we don't care where we are.
      if ( mod.overwriteCrust == true ) 
      {
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
     
      // Here, we do. Multiply by the taper in depth.
      if ( mod.overwriteCrust == false )
      { 
        msh.c11[i] = C      * downTap + msh.c11[i] * upTap;
        msh.c22[i] = A      * downTap + msh.c22[i] * upTap;
        msh.c33[i] = A      * downTap + msh.c33[i] * upTap;
        msh.c12[i] = F      * downTap + msh.c12[i] * upTap;
        msh.c13[i] = F      * downTap + msh.c13[i] * upTap;
        msh.c23[i] = S      * downTap + msh.c23[i] * upTap;
        msh.c44[i] = N      * downTap + msh.c44[i] * upTap;
        msh.c55[i] = L      * downTap + msh.c55[i] * upTap;
        msh.c66[i] = L      * downTap + msh.c66[i] * upTap;
        msh.rho[i] = rhoUse * downTap + msh.rho[i] * upTap;                                                  
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
               double &c56, double &c66, double &rho, char mode )
{
  
  Utilities util;
  Constants con;
 
  // Local trees.
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
  bool found   = false;  
  int  count   = 0;
  int  nodeNum = 0;

  /* *point* stores the current point. *pointClose* stores the closest ever
   * point */
  int  point, pointClose;

  /* These things iterate over the full mesh search */
  int xx0, xx1, xx2, xx3;
  int allIter = 0;

  // Find node closest to point.
  set   = kd_nearest3 ( tree, testX, testY, testZ );
  ind_p = kd_res_item_data ( set );
  point = * ( int * ) ind_p;
  kd_res_free ( set );

  // Store the closest point.
  pointClose = point;
  
  // Find the originally requested col, lon, and rad.
  util.xyz2ColLonRadRad ( origX, origY, origZ, col, lon, rad );        
  
  // Find ColLonRad of original node.
  util.xyz2ColLonRadRad ( msh.xmsh[point], msh.ymsh[point], msh.zmsh[point], 
    colClose, lonClose, radClose );

  // Repeat until enclosing tet is found.
  while ( found == false ) 
  {   

    // Barycentric coordinates.
    double l1, l2, l3, l4;
  
    // Set up connectivity iterator.
    pair < multimap <int, vector <int> > :: iterator , multimap 
      <int, vector <int> > :: iterator > ext;

    // Extract iterator for current point.
    ext = msh.elemOrder.equal_range (point);
    
    /* Mode 'a' for 'allMesh' */
    if ( mode == 'a' )
    {      
      
      // Convert to barycentric coordinates, and test the l values.
      util.convertBary ( origX, origY, origZ,
        msh.xmsh[xx0], msh.xmsh[xx1], 
        msh.xmsh[xx2], msh.xmsh[xx3],
        msh.ymsh[xx0], msh.ymsh[xx1],
        msh.ymsh[xx2], msh.ymsh[xx3],
        msh.zmsh[xx0], msh.zmsh[xx1],
        msh.zmsh[xx2], msh.zmsh[xx3],
        l1, l2, l3, l4 ); 

      if ( l1 >= 0 && l2 >= 0 && l3 >= 0 && l4 >= 0 ) 
      {

        found = true;

        double c11p0 = msh.c11[xx0];
        double c12p0 = msh.c12[xx0];
        double c13p0 = msh.c13[xx0];
        double c14p0 = msh.c14[xx0];
        double c15p0 = msh.c15[xx0];
        double c16p0 = msh.c16[xx0];
        double c22p0 = msh.c22[xx0];
        double c23p0 = msh.c23[xx0];
        double c24p0 = msh.c24[xx0];
        double c25p0 = msh.c25[xx0];
        double c26p0 = msh.c26[xx0];
        double c33p0 = msh.c33[xx0];
        double c34p0 = msh.c34[xx0];
        double c35p0 = msh.c35[xx0];
        double c36p0 = msh.c36[xx0];
        double c44p0 = msh.c44[xx0];
        double c45p0 = msh.c45[xx0];
        double c46p0 = msh.c46[xx0];
        double c55p0 = msh.c55[xx0];
        double c56p0 = msh.c56[xx0];
        double c66p0 = msh.c66[xx0]; 
        double rhop0 = msh.rho[xx0];  
      
        double c11p1 = msh.c11[xx1];
        double c12p1 = msh.c12[xx1];
        double c13p1 = msh.c13[xx1];
        double c14p1 = msh.c14[xx1];
        double c15p1 = msh.c15[xx1];
        double c16p1 = msh.c16[xx1];
        double c22p1 = msh.c22[xx1];
        double c23p1 = msh.c23[xx1];
        double c24p1 = msh.c24[xx1];
        double c25p1 = msh.c25[xx1];
        double c26p1 = msh.c26[xx1];
        double c33p1 = msh.c33[xx1];
        double c34p1 = msh.c34[xx1];
        double c35p1 = msh.c35[xx1];
        double c36p1 = msh.c36[xx1];
        double c44p1 = msh.c44[xx1];
        double c45p1 = msh.c45[xx1];
        double c46p1 = msh.c46[xx1];
        double c55p1 = msh.c55[xx1];
        double c56p1 = msh.c56[xx1];
        double c66p1 = msh.c66[xx1];
        double rhop1 = msh.rho[xx1];     

        double c11p2 = msh.c11[xx2];
        double c12p2 = msh.c12[xx2];
        double c13p2 = msh.c13[xx2];
        double c14p2 = msh.c14[xx2];
        double c15p2 = msh.c15[xx2];
        double c16p2 = msh.c16[xx2];
        double c22p2 = msh.c22[xx2];
        double c23p2 = msh.c23[xx2];
        double c24p2 = msh.c24[xx2];
        double c25p2 = msh.c25[xx2];
        double c26p2 = msh.c26[xx2];
        double c33p2 = msh.c33[xx2];
        double c34p2 = msh.c34[xx2];
        double c35p2 = msh.c35[xx2];
        double c36p2 = msh.c36[xx2];
        double c44p2 = msh.c44[xx2];
        double c45p2 = msh.c45[xx2];
        double c46p2 = msh.c46[xx2];
        double c55p2 = msh.c55[xx2];
        double c56p2 = msh.c56[xx2];
        double c66p2 = msh.c66[xx2];    
        double rhop2 = msh.rho[xx2];     

        double c11p3 = msh.c11[xx3];
        double c12p3 = msh.c12[xx3];
        double c13p3 = msh.c13[xx3];
        double c14p3 = msh.c14[xx3];
        double c15p3 = msh.c15[xx3];
        double c16p3 = msh.c16[xx3];
        double c22p3 = msh.c22[xx3];
        double c23p3 = msh.c23[xx3];
        double c24p3 = msh.c24[xx3];
        double c25p3 = msh.c25[xx3];
        double c26p3 = msh.c26[xx3];
        double c33p3 = msh.c33[xx3];
        double c34p3 = msh.c34[xx3];
        double c35p3 = msh.c35[xx3];
        double c36p3 = msh.c36[xx3];
        double c44p3 = msh.c44[xx3];
        double c45p3 = msh.c45[xx3];
        double c46p3 = msh.c46[xx3];
        double c55p3 = msh.c55[xx3];
        double c56p3 = msh.c56[xx3];
        double c66p3 = msh.c66[xx3];    
        double rhop3 = msh.rho[xx3];   

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

    if ( mode == 'p' )
    {
      // Loop over connecting elements (node indices contained in iterator).
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
          
        // If barycentric coordinates are all >= 0.
        if ( l1 >= 0 && l2 >= 0 && l3 >= 0 && l4 >= 0 ) 
        {
        
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
    }
    
    /* If we haven't found the point on the first try, need to do some magical
     * stuff */
    if ( found == false ) 
    {   
      
      /* First, iterate the count parameter, which keeps track of all the false
       * hits we get. */
      count++;
      
      /* Take the test value, and get the col, lon, and rad of that
       * point */
      util.xyz2ColLonRadRad ( testX, testY, testZ, colPoint, lonPoint, 
        radPoint );            

      /* As an added thing, this checks if we're in some floating underflow
       * off the mesh type thing. It adjusts the originally requested point 
       * to trick the interpolater that the point is inside the mesh 
       * by perturbing its position sligtly). */
      double colOrig, lonOrig, radOrig;      
      util.xyz2ColLonRadDeg ( origX, origY, origZ, colOrig, lonOrig, 
        radOrig );              
      utl.checkMeshEdge     ( colOrig, lonOrig, msh );        
      util.colLonRadDeg2xyz ( colOrig, lonOrig, radOrig, origX, origY, origZ );
      
      /* For col and lon, randomly choose which direction to look. This might
      be able to be switched to a more direction search (i.e. we look in the
      direction which is towards the requested point), but this is certainly
      a bit more robust. For radius, we do a directional search (i.e. search 
      only down if we are higher in radius). This is because radial 
      discretization is much finer than lat/lon. */
      int signC = rand () % 2 ? 1 : -1;  
      int signL = rand () % 2 ? 1 : -1;
      int signR = radClose < rad ? 1 : -1;
      
      /* The theta discretization is controlled by the average edge length of an
      element in the lat/lon direction. These parameters could be optimized 
      significantly. */      
      double dTheta = ( 100 / rad );

      /* Since the size of a lon degree depends so much on the latitude, make
       * this parameter dependent on the latitude. Should change this to some 
       * sort of linear scaling */
      if ( col < 10 )
      {
        dTheta = 170 / rad;
      }

      if ( col < 5 )
      {
        dTheta = 500 / rad;
      }
    
      /* Col and Lon are allowed to vary between their original values (0) and
       * one edge length away. This seems to work. Radius is allowed to vary by
       * 1 km. This also seems to work, although I think it's a bit sketchier. 
       * Could add a parameter to the radius search to make this variable, 
       * dependent on depth. It does work quite well now though. */
      double randC = (rand () % 100) / 100.;
      double randL = (rand () % 100) / 100.;
      double randR = (rand () % 100) / 100.;

      /* Deep mesh is big. This parameter can be changed. Brute force safety 
       * allows a lot of games to be played. */
      if ( rad <= 5371 )
      {
        randR = (rand () % 5000) / 100.;
      }
      
      /* Do *count* iterations of a 'smart' search. Keep looking through the 
       * KDtree with a new, randomly chosen point. */      
      if ( count <= fallBackCount ) 
      {
        double colTest = col + ( signC * randC * dTheta );
        double lonTest = lon + ( signL * randL * dTheta );
        double radTest = rad + ( signR * randR );
        
        /* Create a new testX, Y, and Z, point for the recursive search. */
        util.colLonRadRad2xyz ( colTest, lonTest, radTest, testX, testY, 
          testZ );          
          
        // Extract point from KDTree.
        set   = kd_nearest3 ( tree, testX, testY, testZ );
        ind_p = kd_res_item_data ( set );
        point = * ( int * ) ind_p;  
        
        kd_res_free ( set );
      }

      /* Give the random algorithm *count* times to find the enclosing tet. 
       * This is relatively arbitrary. Performance may increase/decrease by 
       * adjusting this parameter. Switch mode to 'a' and do a fallback
       * complete mesh search. */
      if ( count > fallBackCount )
      {
        mode = 'a';
        xx0  = msh.masterElemConn[allIter+0];
        xx1  = msh.masterElemConn[allIter+1];
        xx2  = msh.masterElemConn[allIter+2];
        xx3  = msh.masterElemConn[allIter+3];

        // Iterate by number of nodes per element.
        allIter = allIter + 4;
      }

      /* If the brute force search fails, just take closest point. This means 
       * that the point is not actually located in the mesh (i.e. floating 
       * underflow has resulted it in being just outside. */
      if ( allIter == (msh.masterElemConn.size()) )
      {

        double colPointDeg = 0;
        double lonPointDeg = 0;
        double radPointDeg = 0;

        util.xyz2ColLonRadDeg ( origX, origY, origZ, colPointDeg, lonPointDeg, 
          radPointDeg ); 
          
        std::cout << "Extracting closest. [col, lon, deg] " << colPointDeg 
          << ' ' << lonPointDeg << ' ' << radPointDeg << ' ' 
          << std::flush << std::endl;

        c11 = msh.c11[pointClose];
        c12 = msh.c12[pointClose];
        c13 = msh.c13[pointClose];
        c14 = msh.c14[pointClose];
        c15 = msh.c15[pointClose];
        c16 = msh.c16[pointClose];
        c22 = msh.c22[pointClose];
        c23 = msh.c23[pointClose];
        c24 = msh.c24[pointClose];
        c25 = msh.c25[pointClose];
        c26 = msh.c26[pointClose];
        c33 = msh.c33[pointClose];
        c34 = msh.c34[pointClose];
        c35 = msh.c35[pointClose];
        c36 = msh.c36[pointClose];
        c44 = msh.c44[pointClose];
        c45 = msh.c45[pointClose];
        c46 = msh.c46[pointClose];
        c55 = msh.c55[pointClose];
        c56 = msh.c56[pointClose];
        c66 = msh.c66[pointClose];    
        rho = msh.rho[pointClose];    
        break;

      }
    }
  }  
   
  // Extract just the index of the closet node if we're looking for it.     
  if ( mode == 'e' ) 
  {
    return msh.node_num_map[point];
  } 
  else 
  { 
    return 0;
  }
  
}  
