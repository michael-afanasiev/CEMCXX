#include <iostream>
#include <assert.h>
#include <ctype.h>

#include "kdtree.h"
#include "classes.hpp"

using namespace std;

void Mesh::getInfo (int exoid) 
{
    
  float dum1;
  char  dum2;
  
  ier = ex_inquire ( exoid, EX_INQ_NODES, &num_nodes, &dum1, &dum2 );
  ier = ex_inquire ( exoid, EX_INQ_ELEM,  &num_elem,  &dum1, &dum2 ); 
  
}

void Mesh::populateParams ( int exoid, Model_file &mod ) 
{
  
  if ( mod.intentions == "INITIZALIZE" ) {
  ier = ex_get_nodal_var ( exoid, 0, 0, num_nodes, c11 );
  }
  
}

void Mesh::populateCoord ( int exoid ) 
{ 
       
  ier = ex_get_coord ( exoid, xmsh, ymsh, zmsh );
  
  if (ier != 0) 
  {
    std::cout << "***Fatal error reading in coordinates. Exiting.\n";
    exit (EXIT_FAILURE);
  }
        
}

void Mesh::allocateMesh ( int &num_nodes ) 
{  
  
  xmsh = new double [num_nodes]();
  ymsh = new double [num_nodes]();
  zmsh = new double [num_nodes]();
  
  c11 = new double [num_nodes]();
  c12 = new double [num_nodes]();
  c13 = new double [num_nodes]();
  c14 = new double [num_nodes]();
  c15 = new double [num_nodes]();
  c16 = new double [num_nodes]();
  c22 = new double [num_nodes]();
  c23 = new double [num_nodes]();
  c24 = new double [num_nodes]();
  c25 = new double [num_nodes]();
  c26 = new double [num_nodes]();
  c33 = new double [num_nodes]();
  c34 = new double [num_nodes]();
  c35 = new double [num_nodes]();
  c36 = new double [num_nodes]();
  c44 = new double [num_nodes]();
  c45 = new double [num_nodes]();
  c46 = new double [num_nodes]();
  c55 = new double [num_nodes]();
  c56 = new double [num_nodes]();
  c66 = new double [num_nodes](); 
  
  cout << num_nodes;
   
}

void Mesh::deallocateMesh ()
{
  
  delete [] xmsh;
  delete [] ymsh;
  delete [] zmsh;

  delete [] c11;
  delete [] c12;
  delete [] c13;
  delete [] c14;
  delete [] c15;
  delete [] c16;
  delete [] c22;
  delete [] c23;
  delete [] c24;
  delete [] c25;
  delete [] c26;
  delete [] c33;
  delete [] c34;
  delete [] c35;
  delete [] c36;
  delete [] c44;
  delete [] c45;
  delete [] c46;
  delete [] c55;
  delete [] c56;
  delete [] c66;
    
}

void Mesh::interpolateModel ( Model_file &mod )
{
  
  Utilities util;
  Constants con;
  double l1;
  double l2;
  double l3;
  double l4;
  
  
  // Create KDtree.
  kdtree *tree = kd_create (3);  
  
  // Create reference array.
  int *dat = new int [mod.num_p];  
  for ( int i=0; i<mod.num_p; i++ ) {
    dat[i] = i;
  }
  
  // Populate KDtree.
  cout << "Creating KDtree.\n";
  for ( int i=0; i<mod.num_p; i++ ) {
    kd_insert3 ( tree, mod.x[i], mod.y[i], mod.z[i], &dat[i] );
    cout << "Creating tree. " << mod.num_p - (i + 1) << " nodes left.\xd";
  }  
  cout << "\n";
  cout << "Done creating tree.\n" << std::flush;
    
  // Loop through model file & assign.
  
  if ( mod.interpolation == "GOURAUD" ) {
    cout << "Performing Gouraud shading.\n";
  }
  
  for ( int i=0; i<num_nodes; i++ ) {      
    
    double mshCol;
    double mshLon;
    double mshRad;
    double modCol;
    double modLon;
    double modRad;
    double xTest;
    double yTest;
    double zTest;
    
    kdres *set  = kd_nearest3 ( tree,  xmsh[i],  ymsh[i], zmsh[i] );    
    void *ind_p = kd_res_item_data ( set );
    int point   = * ( int * ) ind_p;
    kd_res_free ( set );         
    
    double nearX0 = mod.x[point];
    double nearY0 = mod.y[point];
    double nearZ0 = mod.z[point];   
        
    util.xyz2ColLonRadDeg ( mod.x[point], mod.y[point], mod.z[point],
                            modCol,   modLon,   modRad );
    util.xyz2ColLonRadDeg ( xmsh[i], ymsh[i], zmsh[i],
                            mshCol,      mshLon,      mshRad );
    
    int colDir = modCol < mshCol ? 1 : -1;
    int lonDir = modLon < mshLon ? 1 : -1;
    int radDir = modRad < mshRad ? 1 : -1;
    
    int newPointRad = point;
    double testRad  = mshRad;
    while ( newPointRad == point ) {
      
      testRad = testRad + ( radDir * 1 );
      
      util.colLonRadDeg2xyz ( mshCol, mshLon, testRad, 
                              xTest,  yTest,  zTest );
      set         = kd_nearest3 ( tree, xTest, yTest, zTest );
      ind_p       = kd_res_item_data ( set );
      newPointRad = * ( int * ) ind_p;
      kd_res_free ( set );
    
      if ( testRad > con.R_EARTH ) {
        bool surface = true;
        break;
      }
      
    }
    
    double nearX1 = mod.x[newPointRad];
    double nearY1 = mod.y[newPointRad];
    double nearZ1 = mod.z[newPointRad];       
    
    int newPointCol = point;
    double testCol  = mshCol;    
    while ( newPointCol == point ) {
      
      testCol = testCol + ( colDir * 0.1 );
      
      util.colLonRadDeg2xyz ( testCol, mshLon, mshRad, 
                              xTest,  yTest,  zTest );
                              
      set         = kd_nearest3 ( tree, xTest, yTest, zTest );
      ind_p       = kd_res_item_data ( set );
      newPointCol = * ( int * ) ind_p;
      kd_res_free ( set );
      
      if ( mod.x[newPointCol] == nearX1 ) {
        newPointCol = point;
      }
      
    }
    
    double nearX2 = mod.x[newPointCol];
    double nearY2 = mod.y[newPointCol];
    double nearZ2 = mod.z[newPointCol];       
    
    int newPointLon = point;
    double testLon  = mshLon;    
    while ( newPointLon == point ) {
      
      testLon = testLon + ( lonDir * 0.1 );
    
      util.colLonRadDeg2xyz ( mshCol, testLon, mshRad, 
                              xTest,  yTest,  zTest );
                              
      set         = kd_nearest3 ( tree, xTest, yTest, zTest );
      ind_p       = kd_res_item_data ( set );
      newPointLon = * ( int * ) ind_p;
      kd_res_free ( set );
      
      if ( mod.x[newPointLon] == nearX2 || mod.x[newPointLon] == nearX1 ) {
        newPointLon = point;
      }
      
    }
    
    double nearX3 = mod.x[newPointLon];
    double nearY3 = mod.y[newPointLon];
    double nearZ3 = mod.z[newPointLon];   
        
    cout << "Requested point: " << xmsh[i] << " " << ymsh[i] << " " << 
      zmsh[i] << "\n";
    cout << "Found point0: " << nearX0 << " " << nearY0 << " " << 
      nearZ0 << "\n";
    cout << "Found point1: " << nearX1 << " " << nearY1 << " " << 
      nearZ1 << "\n";
    cout << "Found point2: " << nearX2 << " " << nearY2 << " " << 
      nearZ2 << "\n";
    cout << "Found point3: " << nearX3 << " " << nearY3 << " " << 
      nearZ3 << "\n";
    
    cin.get ();
    
                      
    util.convertBary ( mod.x[i], mod.y[i], mod.z[i],
                       nearX0, nearX1, nearX2, nearX3,
                       nearY0, nearY1, nearY2, nearY3,
                       nearZ0, nearZ1, nearZ2, nearZ3,
                       l1, l2, l3, l4 );
                       
   // cout << "Barycentric: " << l1 << " " << l2 << " " << l3 << " " << l4 << "\n";
   // cin.get ();
                       
                       
                                         
    cout << "Assinging model. " << mod.x.size() - (i + 1) << " nodes left.\xd";
  }
  cout << "\n";

}

void Mesh::reNormalize ( Model_file &mod ) 
{
  
  Constants con;
  
  if ( mod.dimensions == "NORMALIZED" ) {
    for ( int i=0; i<num_nodes; i++ ) {
      xmsh[i] = xmsh[i] * con.R_EARTH;
      ymsh[i] = ymsh[i] * con.R_EARTH;
      zmsh[i] = zmsh[i] * con.R_EARTH;
    }
  }
    
}