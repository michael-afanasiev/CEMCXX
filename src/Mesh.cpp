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

void Mesh::getInfo ( int in_exoid ) 
{
    
  exoid = in_exoid;  
  
  float dum1;
  char  dum2;
  
  ier = ex_inquire ( exoid, EX_INQ_NODES,    &num_nodes,    &dum1, &dum2 );
  ier = ex_inquire ( exoid, EX_INQ_ELEM,     &num_elem,     &dum1, &dum2 ); 
  ier = ex_inquire ( exoid, EX_INQ_ELEM_BLK, &num_elem_blk, &dum1, &dum2 );
  
  allocateMesh   ();  
  populateParams ();
  populateCoord  ();
  
}

void Mesh::populateParams ( ) 
{

  ier = ex_get_nodal_var ( exoid, 1, 1,  num_nodes, c11 );
  ier = ex_get_nodal_var ( exoid, 1, 2,  num_nodes, c12 );
  ier = ex_get_nodal_var ( exoid, 1, 3,  num_nodes, c13 );
  ier = ex_get_nodal_var ( exoid, 1, 4,  num_nodes, c14 );
  ier = ex_get_nodal_var ( exoid, 1, 5,  num_nodes, c15 );
  ier = ex_get_nodal_var ( exoid, 1, 6,  num_nodes, c16 );
  ier = ex_get_nodal_var ( exoid, 1, 7,  num_nodes, c22 );
  ier = ex_get_nodal_var ( exoid, 1, 8,  num_nodes, c23 );
  ier = ex_get_nodal_var ( exoid, 1, 9,  num_nodes, c24 );
  ier = ex_get_nodal_var ( exoid, 1, 10, num_nodes, c25 );
  ier = ex_get_nodal_var ( exoid, 1, 11, num_nodes, c26 );
  ier = ex_get_nodal_var ( exoid, 1, 12, num_nodes, c33 );
  ier = ex_get_nodal_var ( exoid, 1, 13, num_nodes, c34 );
  ier = ex_get_nodal_var ( exoid, 1, 14, num_nodes, c35 );
  ier = ex_get_nodal_var ( exoid, 1, 15, num_nodes, c36 );
  ier = ex_get_nodal_var ( exoid, 1, 16, num_nodes, c44 );
  ier = ex_get_nodal_var ( exoid, 1, 17, num_nodes, c45 );
  ier = ex_get_nodal_var ( exoid, 1, 18, num_nodes, c46 );
  ier = ex_get_nodal_var ( exoid, 1, 19, num_nodes, c55 );
  ier = ex_get_nodal_var ( exoid, 1, 20, num_nodes, c56 );
  ier = ex_get_nodal_var ( exoid, 1, 21, num_nodes, c66 );
  ier = ex_get_nodal_var ( exoid, 1, 22, num_nodes, rho );
  
  if ( ier != 0 ) {
    cout << "Error reading in mesh variables.\n";
    exit (EXIT_FAILURE);
  }
  
}

void Mesh::populateCoord ( ) 
{ 
       
  ier = ex_get_coord ( exoid, xmsh, ymsh, zmsh );
  
  if (ier != 0) {
    std::cout << "***Fatal error reading in coordinates. Exiting.\n";
    exit (EXIT_FAILURE);
  }
}    

void Mesh::allocateMesh ( ) 
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
  rho = new double [num_nodes]();
   
}

void Mesh::getConnectivity ( int exoid )
{
  
  int *ids = new int [num_elem_blk];
  int ier  = ex_get_elem_blk_ids ( exoid, ids );
  
  int *elemConn = new int [num_elem_blk*num_elem*num_node_per_elem];  
  ier           = ex_get_elem_conn ( exoid, ids[0], elemConn );
  
  vector<int> node;
  node.reserve ( num_node_per_elem );

  cout << "Building connectivity array.\n";  
  for ( int i=0; i<num_elem_blk*num_elem*num_node_per_elem; i++ ) {
    node.push_back ( elemConn[i] - 1 );
    
    if ( (i+1) % num_node_per_elem == 0 ) {            
      elemOrder.insert ( pair <int, vector <int> > (node[0], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[1], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[2], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[3], node) );                               
      
      node.clear ();
    }    
  }
      
  delete [] elemConn;
  
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
  
delete [] rho;
    
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