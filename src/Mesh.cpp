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

void Mesh::getInfo ( int in_exoid, char mode ) 
{
    
  exoid = in_exoid;  
  
  float dum1;
  char  dum2;
  
  ier = ex_inquire ( exoid, EX_INQ_NODES,    &num_nodes,    &dum1, &dum2 );
  ier = ex_inquire ( exoid, EX_INQ_ELEM,     &num_elem,     &dum1, &dum2 ); 
  ier = ex_inquire ( exoid, EX_INQ_ELEM_BLK, &num_elem_blk, &dum1, &dum2 );
  
  allocateMesh   ();  
  populateCoord  ();
  
  if ( mode == 'p' ) {
    populateParams ();
  }
}

void Mesh::populateParams ( ) 
{
  
  int numMshVar = 27;
  
  char *var_names[numMshVar];
  for ( int i=0; i<numMshVar; i++ )
    var_names[i] = (char *) calloc ( (MAX_STR_LENGTH+1), sizeof(char) );
  
  ier = ex_get_var_names ( exoid, "n", 27, var_names );
    
  int c11i, c12i, c13i, c14i, c15i, c16i, c22i, c23i, c24i, c25i, c26i, c33i;
  int c34i, c35i, c36i, c44i, c45i, c46i, c55i, c56i, c66i, rhoi, Q__i, elvi;
  int du1i, du2i, du3i;
  
  for ( int i=0; i<numMshVar; i++ ) {
    
    if (( strcmp ( var_names[i], "c11" ) ) == 0 )
      c11i = i + 1;
    if (( strcmp ( var_names[i], "c12" ) ) == 0 )
      c12i = i + 1;
    if (( strcmp ( var_names[i], "c13" ) ) == 0 )
      c13i = i + 1;
    if (( strcmp ( var_names[i], "c14" ) ) == 0 )
      c14i = i + 1;
    if (( strcmp ( var_names[i], "c15" ) ) == 0 )
      c15i = i + 1;
    if (( strcmp ( var_names[i], "c16" ) ) == 0 )
      c16i = i + 1;
    if (( strcmp ( var_names[i], "c22" ) ) == 0 )
      c22i = i + 1;
    if (( strcmp ( var_names[i], "c23" ) ) == 0 )
      c23i = i + 1;
    if (( strcmp ( var_names[i], "c24" ) ) == 0 )
      c24i = i + 1;
    if (( strcmp ( var_names[i], "c25" ) ) == 0 )
      c25i = i + 1;
    if (( strcmp ( var_names[i], "c26" ) ) == 0 )
      c26i = i + 1;
    if (( strcmp ( var_names[i], "c33" ) ) == 0 )
      c33i = i + 1;
    if (( strcmp ( var_names[i], "c34" ) ) == 0 )
      c34i = i + 1;
    if (( strcmp ( var_names[i], "c35" ) ) == 0 )
      c35i = i + 1;
    if (( strcmp ( var_names[i], "c36" ) ) == 0 )
      c36i = i + 1;
    if (( strcmp ( var_names[i], "c44" ) ) == 0 )
      c44i = i + 1;    
    if (( strcmp ( var_names[i], "c45" ) ) == 0 )
      c45i = i + 1;
    if (( strcmp ( var_names[i], "c46" ) ) == 0 )
      c46i = i + 1;
    if (( strcmp ( var_names[i], "c55" ) ) == 0 )
      c55i = i + 1;
    if (( strcmp ( var_names[i], "c56" ) ) == 0 )
      c56i = i + 1;
    if (( strcmp ( var_names[i], "c66" ) ) == 0 )
      c66i = i + 1;
    if (( strcmp ( var_names[i], "rho" ) ) == 0 )
      rhoi = i + 1;
    if (( strcmp ( var_names[i], "Q__" ) ) == 0 )
      Q__i = i + 1;
    if (( strcmp ( var_names[i], "elv" ) ) == 0 )
      elvi = i + 1;
    if (( strcmp ( var_names[i], "du1" ) ) == 0 )
      du1i = i + 1;
    if (( strcmp ( var_names[i], "du2" ) ) == 0 )
      du2i = i + 1;    
    if (( strcmp ( var_names[i], "du3" ) ) == 0 )
      du3i = i + 1;
    
  }  
  
  ier = ex_get_nodal_var ( exoid, 1, c11i, num_nodes, c11 );
  ier = ex_get_nodal_var ( exoid, 1, c12i, num_nodes, c12 );
  ier = ex_get_nodal_var ( exoid, 1, c13i, num_nodes, c13 );
  ier = ex_get_nodal_var ( exoid, 1, c14i, num_nodes, c14 );
  ier = ex_get_nodal_var ( exoid, 1, c15i, num_nodes, c15 );
  ier = ex_get_nodal_var ( exoid, 1, c16i, num_nodes, c16 );
  ier = ex_get_nodal_var ( exoid, 1, c22i, num_nodes, c22 );
  ier = ex_get_nodal_var ( exoid, 1, c23i, num_nodes, c23 );
  ier = ex_get_nodal_var ( exoid, 1, c24i, num_nodes, c24 );
  ier = ex_get_nodal_var ( exoid, 1, c25i, num_nodes, c25 );
  ier = ex_get_nodal_var ( exoid, 1, c26i, num_nodes, c26 );
  ier = ex_get_nodal_var ( exoid, 1, c33i, num_nodes, c33 );
  ier = ex_get_nodal_var ( exoid, 1, c34i, num_nodes, c34 );
  ier = ex_get_nodal_var ( exoid, 1, c35i, num_nodes, c35 );
  ier = ex_get_nodal_var ( exoid, 1, c36i, num_nodes, c36 );
  ier = ex_get_nodal_var ( exoid, 1, c44i, num_nodes, c44 );
  ier = ex_get_nodal_var ( exoid, 1, c45i, num_nodes, c45 );
  ier = ex_get_nodal_var ( exoid, 1, c46i, num_nodes, c46 );
  ier = ex_get_nodal_var ( exoid, 1, c55i, num_nodes, c55 );
  ier = ex_get_nodal_var ( exoid, 1, c56i, num_nodes, c56 );
  ier = ex_get_nodal_var ( exoid, 1, c66i, num_nodes, c66 );
  ier = ex_get_nodal_var ( exoid, 1, rhoi, num_nodes, rho );
  ier = ex_get_nodal_var ( exoid, 1, Q__i, num_nodes, Q__ );
  ier = ex_get_nodal_var ( exoid, 1, elvi, num_nodes, elv );
  ier = ex_get_nodal_var ( exoid, 1, du1i, num_nodes, du1 );
  ier = ex_get_nodal_var ( exoid, 1, du2i, num_nodes, du2 );
  ier = ex_get_nodal_var ( exoid, 1, du3i, num_nodes, du3 );
    
  /* FIXME There is a bug in the getting of the rho variable. It looks like all
  variables need to be filled when writing an exodus file, otherwise the last 
  variable will be written at the end. So, elv, dum1, dum2, and dum3 need to be
  written to make this work properly. */  
  
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
  Q__ = new double [num_nodes]();
  elv = new double [num_nodes]();
  du1 = new double [num_nodes]();
  du2 = new double [num_nodes]();
  du3 = new double [num_nodes]();
     
}

void Mesh::getConnectivity ( int exoid )
{
        
  int *ids = new int [num_elem_blk];
  int ier  = ex_get_elem_blk_ids ( exoid, ids );
  
  vector <int> masterElemConn;
  for ( int i=0; i!=num_elem_blk; i++ ) {
    
    ier = ex_get_elem_block ( exoid, ids[i], elem_type, &num_elem_in_blk, 
      &num_nodes_in_elem, &num_attr );
  
    int *elemConn = new int [num_elem_in_blk*num_node_per_elem];  
    ier           = ex_get_elem_conn ( exoid, ids[i], elemConn );
    
    for ( int j=0; j!= num_elem_in_blk*num_node_per_elem; j++ ) {
      masterElemConn.push_back (elemConn[j]);
    }
    
    delete [] elemConn;
    
  }
    
  vector <int> node;
  node.reserve ( num_node_per_elem );
  
  cout << "Building connectivity array.\n";  
  for ( int i=0; i<num_elem*num_node_per_elem; i++ ) {
    node.push_back ( masterElemConn[i] - 1 );
    
    if ( (i+1) % num_node_per_elem == 0 ) {            
      elemOrder.insert ( pair <int, vector <int> > (node[0], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[1], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[2], node) );
      elemOrder.insert ( pair <int, vector <int> > (node[3], node) );  
            
      node.clear ();
    }    
  }   
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