#include <iostream>
#include <assert.h>
#include <ctype.h>
#include "kdtree.h"
#include "classes.h"
using namespace std;

void Mesh::getInfo (int exoid) 
{
    
  float dum1;
  char  dum2;
  
  ier = ex_inquire ( exoid, EX_INQ_NODES, &num_nodes, &dum1, &dum2 );
  ier = ex_inquire ( exoid, EX_INQ_ELEM,  &num_elem,  &dum1, &dum2 ); 
  
}

void Mesh::populateParams ( int exoid ) 
{
  
  ier = ex_get_nodal_var ( exoid, 0, 0, num_nodes, c11 );
  
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

  // Create KDtree.
  kdtree *tree = kd_create (3);  
  
  // Create reference array.
  int *dat = new int [num_nodes];  
  for ( int i=0; i<num_nodes; i++ ) {
    dat[i] = i;
  }
  
  // Populate KDtree.
  cout << "Creating KDtree.\n";
  for ( int i=0; i<num_nodes; i++ ) {
    kd_insert3 ( tree, xmsh[i], ymsh[i], zmsh[i], &dat[i] );
    cout << "Creating tree. " << num_nodes - (i + 1) << " nodes left.\xd";
  }  
  cout << "\n";
  cout << "Done creating tree.\n" << std::flush;
  
  // Loop through model file & assign.
  
  if ( mod.interpolation == "GOURAUD" ) {
    cout << "Performing Gouraud shading.\n";
  }
  
  for ( int i=0; i<mod.x.size(); i++ ) {      
    kdres  *set = kd_nearest3  ( tree,  mod.x[i],  mod.y[i], mod.z[i] );
    void *ind_p = kd_res_item_data ( set );
    int     ind = * (int *) ind_p;

    cout << "Assinging model. " << mod.x.size() - (i + 1) << " nodes left.\xd";
    kd_res_free ( set );        
  }
  cout << "\n";

}