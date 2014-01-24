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

  int *ids = new int [num_elem_blk];
  int ier  = ex_get_elem_blk_ids ( exoid, ids );
  
  int *elemConn = new int [num_elem_blk*num_elem*num_nodes];  
  ier           = ex_get_elem_conn ( exoid, ids[0], elemConn );
  
  elemOrder.resize ( num_nodes );
  
  vector<int> node;
  node.reserve ( 4 );
  
  for ( int i=0; i<num_elem_blk*num_elem; i++ ) {
    node.push_back ( elemConn[i] );
    if ( (i+1) % 4 == 0 ) {
      elemOrder[node[0]-1].push_back (node[0]);
      elemOrder[node[0]-1].push_back (node[1]);
      elemOrder[node[0]-1].push_back (node[2]);
      elemOrder[node[0]-1].push_back (node[3]);
                       
      elemOrder[node[1]-1].push_back (node[0]);
      elemOrder[node[1]-1].push_back (node[1]);
      elemOrder[node[1]-1].push_back (node[2]);
      elemOrder[node[1]-1].push_back (node[3]);
                       
      elemOrder[node[2]-1].push_back (node[0]);
      elemOrder[node[2]-1].push_back (node[1]);
      elemOrder[node[2]-1].push_back (node[2]);
      elemOrder[node[2]-1].push_back (node[3]);
                       
      elemOrder[node[3]-1].push_back (node[0]);
      elemOrder[node[3]-1].push_back (node[1]);
      elemOrder[node[3]-1].push_back (node[2]);
      elemOrder[node[3]-1].push_back (node[3]);
                        
      
      cout << "elemOrder: " << elemOrder[node[0]][0] << "\n";
      cout << "elemOrder: " << elemOrder[node[0]][1] << "\n";
      cout << "elemOrder: " << elemOrder[node[0]][2] << "\n";
      cout << "elemOrder: " << elemOrder[node[0]][3] << "\n";
      cout << "elemOrder: " << elemOrder[node[0]][4] << "\n";
      cout << "elemOrder: " << elemOrder[node[0]][5] << "\n";
      cout << "elemOrder: " << elemOrder[node[0]][6] << "\n";
      cout << "elemOrder: " << node[0]-1 << "\n";
      
      cout << "node: " << node[0] << "\n";
      cout << "node: " << node[1] << "\n";
      cout << "node: " << node[2] << "\n";
      cout << "node: " << node[3] << "\n";
      
      node.clear ();
      cin.get ();
    }
  }
    
  // list <int> table [num_nodes];
  // list <int>::iterator it;
  // for ( int i=0; i<num_elem_blk*num_elem; i++ ) {
  //   table[0].push_back (elemOrder[i][0]);
  //   table[0].push_back (elemOrder[i][1]);
  //   table[0].push_back (elemOrder[i][2]);
  //   table[0].push_back (elemOrder[i][3]); 
  //         
  //   table[0].push_back (elemOrder[i][0]);
  //   table[0].push_back (elemOrder[i][1]);
  //   table[0].push_back (elemOrder[i][2]);
  //   table[0].push_back (elemOrder[i][3]);  
  //         
  //   table[0].push_back (elemOrder[i][0]);
  //   table[0].push_back (elemOrder[i][1]);
  //   table[0].push_back (elemOrder[i][2]);
  //   table[0].push_back (elemOrder[i][3]);  
  //         
  //   table[0].push_back (elemOrder[i][0]);
  //   table[0].push_back (elemOrder[i][1]);
  //   table[0].push_back (elemOrder[i][2]);
  //   table[0].push_back (elemOrder[i][3]);    
  //   
  //   // for ( it=table[0].begin(); it!=table[0].end(); it++ )
  //   //   cout << " " << setw(5) << *it;
  //   cin.get ();
  // }


  // for ( int i=0; i<num_elem_blk*num_elem; i++ ) {
  //   sort ( elemOrder[i].begin(), elemOrder[i].end() );
  // }
  // stable_sort ( elemOrder.begin(), elemOrder.end() );
  // 
  // int *nodePointer = new int [num_elem*num_nodes];
  // int j            = 1;
  // nodePointer[0]   = 0;
  // for ( int i=0; i<num_elem_blk*num_elem-1; i++ ) {
  //   // if ( elemOrder[i+1][0] > elemOrder[i][0] ) {
  //     cout << "elemOrder: " << elemOrder[i][0] << "\n";
  //     cout << "elemOrder: " << elemOrder[i][1] << "\n";
  //     cout << "elemOrder: " << elemOrder[i][2] << "\n";
  //     cout << "elemOrder: " << elemOrder[i][3] << "\n";
  //     cin.get ();
  //     nodePointer[j] = i + 1;
  //     j++;
  //   // }
  // }
  
  delete [] elemConn;

  // cout << "elemOrder: " << elemOrder[3][0] << "\n";
  // cout << "elemOrder: " << elemOrder[3][1] << "\n";
  // cout << "elemOrder: " << elemOrder[3][2] << "\n";
  // cout << "elemOrder: " << elemOrder[3][3] << "\n";
  // 
    
  // cout << "elemOrder: " << elemOrder[3][0] << "\n";
  // cout << "elemOrder: " << elemOrder[3][1] << "\n";
  // cout << "elemOrder: " << elemOrder[3][2] << "\n";
  // cout << "elemOrder: " << elemOrder[3][3] << "\n";
  // cout << "EXIT\n";
  // exit (EXIT_SUCCESS);
  
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