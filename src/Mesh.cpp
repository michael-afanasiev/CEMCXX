#include <iostream>
#include <algorithm>
#include <vector>
#include <list>
#include <iomanip>
#include <assert.h>
#include <ctype.h>

#include "kdtree.h"
#include "classes.hpp"

using namespace std;

void Mesh::getInfo (int exoid) 
{
    
  float dum1;
  char  dum2;
  
  ier = ex_inquire ( exoid, EX_INQ_NODES,    &num_nodes,    &dum1, &dum2 );
  ier = ex_inquire ( exoid, EX_INQ_ELEM,     &num_elem,     &dum1, &dum2 ); 
  ier = ex_inquire ( exoid, EX_INQ_ELEM_BLK, &num_elem_blk, &dum1, &dum2);
  
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

void Mesh::getConnectivity ( int exoid )
{
  
  cout << "Building connectivity array.\n";

  int *ids = new int [num_elem_blk];
  int ier  = ex_get_elem_blk_ids ( exoid, ids );
  
  int *elemConn = new int [num_elem_blk*num_elem*4];  
  ier           = ex_get_elem_conn ( exoid, ids[0], elemConn );
  
  vector<int> node;
  node.reserve ( 4 );
  
  for ( int i=0; i<num_elem_blk*num_elem*4; i++ ) {
    node.push_back ( elemConn[i] );
    
    if ( (i+1) % 4 == 0 ) {
      elemOrder.insert ( pair <int, vector <int> > (node[0], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[1], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[2], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[3], node) );                               
      
      node.clear ();
    }
  }
      
  delete [] elemConn;
  
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