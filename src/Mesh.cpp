#include <algorithm>
#include <assert.h>
#include <ctype.h>
#include <cstring>

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
  getMinMaxRad   ();
  
  if ( mode == 'p' ) {
    populateParams ();
  }
}

void Mesh::createKDTreeUnpacked ( )
{
 
  cout << "Creating KDTree ( mesh ).\n";
  tree     = kd_create (3);
  KDdat    = new int [num_nodes];
  for ( int i=0; i<num_nodes; i++ ) 
  {
    KDdat[i] = i;
    kd_insert3 ( tree, xmsh[i], ymsh[i], zmsh[i], &KDdat[i] );
  }
  
}

void Mesh::getMinMaxRad ( )
{
  
  Utilities util;
  Constants con;
  
  double col, lon, rad;
  
  cout << "Determining col/lon/rad box.\n";
  
  for ( int i=0; i<num_nodes; i++ ) {
    
    util.xyz2ColLonRadDeg ( xmsh[i], ymsh[i], zmsh[i], col, lon, rad );    
            
    if ( (xmsh[i] != 0) || (ymsh[i] != 0) ) 
    {
      if ( col < colMin )
        colMin = col;
      if ( col > colMax )
        colMax = col;
      if ( lon < lonMin )
        lonMin = lon;
      if ( lon > lonMax )
        lonMax = lon;
      if ( rad < radMin )
        radMin = rad;
      if ( rad > radMax )
        radMax = rad;          
    }
    
    if ( (xmsh[i] == 0) && (ymsh[i] == 0) && (zmsh[i] > 0) )
      colMin = 0;
    if ( (xmsh[i] == 0) && (ymsh[i] == 0) && (zmsh[i] < 0) )
      colMax = 180;
    
  }
          
  if ( (lonMax == 180.) && (lonMin < 0.) ) 
  {    
    lonMin = -90;
    lonMax = -180;
  }
  
  if ( lonMax < lonMin )
    swap ( lonMax, lonMin );
        
  colMin = colMin * con.PI / con.o80;
  colMax = colMax * con.PI / con.o80;
  lonMin = lonMin * con.PI / con.o80;
  lonMax = lonMax * con.PI / con.o80;   
  
  if ( radMax <= (con.innerCoreRad + 1) )
    regString = "innerCore";
  if ( (radMax <= (con.outerCoreRad + 1)) && (radMin >= (con.innerCoreRad - 1)) )
    regString = "outerCore";
  if ( (radMax <= (con.R670 + 1)) && (radMin >= (con.outerCoreRad - 1)) )
    regString = "lowerMantle";
  if ( (radMax <= (con.R400 + 1)) && (radMin >= (con.R670 - 1)) )
    regString = "transitionZone";
  if ( radMin >= (con.R400 - 1) )
    regString = "upperMantle";
       
}

void Mesh::populateParams ( ) 
{
  
  int numMshVar = 28;
  
  char *var_names[numMshVar];
  for ( int i=0; i<numMshVar; i++ )
    var_names[i] = (char *) calloc ( (MAX_STR_LENGTH+1), sizeof(char) );
  
  ier = ex_get_var_names ( exoid, "n", numMshVar, var_names );
    
  int c11i = 0; 
  int c12i = 0; 
  int c13i = 0; 
  int c14i = 0; 
  int c15i = 0; 
  int c16i = 0; 
  int c22i = 0; 
  int c23i = 0; 
  int c24i = 0; 
  int c25i = 0; 
  int c26i = 0; 
  int c33i = 0;
  int c34i = 0; 
  int c35i = 0; 
  int c36i = 0; 
  int c44i = 0; 
  int c45i = 0; 
  int c46i = 0; 
  int c55i = 0; 
  int c56i = 0; 
  int c66i = 0; 
  int rhoi = 0; 
  int Q__i = 0; 
  int elvi = 0;
  int du1i = 0; 
  int du2i = 0; 
  int du3i = 0; 
  int sizi = 0;
  
  for ( int i=0; i<numMshVar; i++ ) 
  {  
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
    if (( strcmp ( var_names[i], "siz" ) ) == 0 )
      sizi = i + 1;    
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
  ier = ex_get_nodal_var ( exoid, 1, sizi, num_nodes, siz );
  
  if ( ier != 0 ) 
  {
    cout << "***";
    cout << "Error reading in mesh variables.";
    cout << "Don't worry, this is ok if you're refining." << endl;
  }    
  
}

void Mesh::populateCoord ( ) 
{ 
       
  ier = ex_get_coord ( exoid, xmsh, ymsh, zmsh );    
  
  if (ier != 0) 
  {
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
  siz = new double [num_nodes]();
     
}

void Mesh::getNodeNumMap   ( int exoid )
{
    
  node_num_map = new int [num_nodes];
  ier = ex_get_node_num_map ( exoid, node_num_map );
  
}

void Mesh::getElemNumMap ( int exoid )
{
  
  elem_num_map = new int [num_elem];
  ier = ex_get_elem_num_map ( exoid, elem_num_map );
  
}

void Mesh::getElementConnectivity ( int exoid, Model_file &mod )
{
  
  Utilities utl;
  
  int *ids = new int [num_elem_blk];
  int ier  = ex_get_elem_blk_ids ( exoid, ids );
  
  if ( ier != 0 ) 
  {
    cout << "Error in getting elem block ids" << endl;
    exit ( EXIT_FAILURE );
  }
    
  for ( int i=0; i!=num_elem_blk; i++ ) 
  {
    
    ier = ex_get_elem_block ( exoid, ids[i], elem_type, &num_elem_in_blk, 
      &num_nodes_in_elem, &num_attr );

    if ( ier != 0 ) 
    {
      cout << "Error in getting elem block data" << endl;
      exit ( EXIT_FAILURE );
    }
  
    std::vector <int> tmpElemConn;
    std::vector <int> tmpNewElemConn;
    
    tmpElemConn.reserve    ( num_elem_in_blk*num_node_per_elem );
    tmpNewElemConn.reserve ( num_elem_in_blk*num_node_per_elem );
      
    int *elemConn = new int [num_elem_in_blk*num_node_per_elem];  
    ier           = ex_get_elem_conn ( exoid, ids[i], elemConn );
    
    if ( ier != 0 ) 
    {
      cout << "Error in getting elem connectivity data" << endl;
      exit ( EXIT_FAILURE );
    }
        
    int outCount = 0;
    int inCount  = 0;    
    for ( int j=0; j<num_elem_in_blk*num_node_per_elem; j++ ) 
    {

      // Define local variables.
      double mshColRot, mshLonRot, mshRadRot; // Sph. Coord. in rotated domain.
      double xRot,      yRot,      zRot;      // Cart. Coord in rotated domain.
      
      int ind = elemConn[j]-1;
      utl.rotateBackward ( xmsh[ind], ymsh[ind], zmsh[ind], xRot, yRot, zRot, 
        mod );
      utl.xyz2ColLonRadDeg ( xRot, yRot, zRot, mshColRot, mshLonRot, 
        mshRadRot );
      
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
        if ( inCount == num_node_per_elem )
        {
          tmpNewElemConn.push_back ( elemConn[j-3] );
          tmpNewElemConn.push_back ( elemConn[j-2] );
          tmpNewElemConn.push_back ( elemConn[j-1] );
          tmpNewElemConn.push_back ( elemConn[j]   );
        }
      }            
      else
      {
        inCount++;
        if ( inCount == num_node_per_elem )
        {
          tmpElemConn.push_back ( elemConn[j-3] );
          tmpElemConn.push_back ( elemConn[j-2] );
          tmpElemConn.push_back ( elemConn[j-1] );
          tmpElemConn.push_back ( elemConn[j]   );                              
        }
      
      }
      
      if ( (outCount % num_node_per_elem) == 0 )
      {
        inCount = 0;
      }            
          
    }
    
    if ( tmpElemConn.size() != 0 )
    {
      refineElemConn.push_back ( tmpElemConn );
    }
    if ( tmpNewElemConn.size() != 0 )
    {
      refineElemConn.push_back ( tmpNewElemConn );
    }
    
    tmpElemConn.clear();      
    tmpNewElemConn.clear();
    delete [] elemConn;
  }  
}

void Mesh::getConnectivity ( int exoid )
{
        
  int *ids = new int [num_elem_blk];
  int ier  = ex_get_elem_blk_ids ( exoid, ids );      
  
  if ( ier != 0 ) 
  {
    cout << "Error in getting elem block ids" << endl;
    exit ( EXIT_FAILURE );
  }
  
  vector <int> masterElemConn;
  for ( int i=0; i!=num_elem_blk; i++ ) {
    
    ier = ex_get_elem_block ( exoid, ids[i], elem_type, &num_elem_in_blk, 
      &num_nodes_in_elem, &num_attr );
      
    if ( ier != 0 ) 
    {
      cout << "Error in getting elem block data" << endl;
      exit ( EXIT_FAILURE );
    }
  
    int *elemConn = new int [num_elem_in_blk*num_node_per_elem];  
    ier           = ex_get_elem_conn ( exoid, ids[i], elemConn );
    
    if ( ier != 0 ) 
    {
      cout << "Error in getting elem connectivity data" << endl;
      exit ( EXIT_FAILURE );
    }
    
    
    for ( int j=0; j!= num_elem_in_blk*num_node_per_elem; j++ ) {
      masterElemConn.push_back (elemConn[j]);
    }
    
    delete [] elemConn;
    
  }
    
  vector <int> node;
  //vector < vector <int> > elemOrderVec;
  node.reserve ( num_node_per_elem );

  elemOrderVec.reserve ( num_elem*num_node_per_elem );
  for ( int isize=0; isize!=num_elem*num_node_per_elem-1; isize++ )
  { 
//    std::cout << isize << ' ' << num_elem*num_node_per_elem << std::flush << std::endl;
    elemOrderVec[isize].reserve(num_node_per_elem);
  }
  
  cout << "Building connectivity array. Done allocating" << std::flush << std::endl; 
  
  for ( int i=0; i<num_elem*num_node_per_elem; i++ ) 
  {
    node.push_back ( masterElemConn[i] - 1 );
    
    if ( (i+1) % num_node_per_elem == 0 ) 
    {            
      elemOrderVec[node[0]].push_back ( node[0] );
      elemOrderVec[node[0]].push_back ( node[1] );
      elemOrderVec[node[0]].push_back ( node[2] );
      elemOrderVec[node[0]].push_back ( node[3] );
      elemOrderVec[node[1]].push_back ( node[0] );
      elemOrderVec[node[1]].push_back ( node[1] );
      elemOrderVec[node[1]].push_back ( node[2] );
      elemOrderVec[node[1]].push_back ( node[3] );
      elemOrderVec[node[2]].push_back ( node[0] );
      elemOrderVec[node[2]].push_back ( node[1] );
      elemOrderVec[node[2]].push_back ( node[2] );
      elemOrderVec[node[2]].push_back ( node[3] );
      elemOrderVec[node[3]].push_back ( node[0] );  
      elemOrderVec[node[3]].push_back ( node[1] );  
      elemOrderVec[node[3]].push_back ( node[2] );  
      elemOrderVec[node[3]].push_back ( node[3] );  
            
      node.clear ();
    }    
  }     
}


void Mesh::deallocateMesh ( Model_file &mod )
{
  
  cout << "Debug 1" << endl;
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
  delete [] Q__;
  delete [] elv;
  delete [] du1;
  delete [] du2;
  delete [] du3;
  delete [] siz;
 
  cout << "Debug 2" << endl;

  if ( mod.intentions == "EXTRACT" )
  {    
    delete [] KDdat;
    kd_free ( tree );  
  }
 
  cout << "Debug 3" << endl;
}
