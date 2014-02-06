#include <stdio.h>
#include <cstdio>
#include <string.h>
#include "classes.hpp"

using namespace std;

void Exodus_file::openFile ( string fname ) 
{
    
  cout << "Opening exodus file: " << fname << "\n";
        
  idexo = ex_open ( fname.c_str(), EX_WRITE, &comp_ws, &io_ws, &vers ); 
  
  if (idexo < 0) {
    std::cout << "***Fatal error opening exodus file. Exiting.\n";  
    exit (EXIT_FAILURE);  
  } 
  
} 

void Exodus_file::closeFile () 
{
  
  ier = ex_close ( idexo );
  
  if (ier == 0) {
    std::cout << "File closed succesfully. \n";
  } else {
    cout << "***Fatal error closing exodus file. Exiting\n";
    exit (EXIT_FAILURE);  
  }
  
}

void Exodus_file::merge ( Model_file &mod ) 
{
  
  
  
}

void Exodus_file::writeParams ( Mesh &msh )
{

  char *cstr = new char [MAX_LINE_LENGTH];
  const char *varnames[22];

  int ndim;
  int nump;
  int numel;
  int numelblk;
  int numnps;
  int numess;
  
  varnames [0]  = "c11";
  varnames [1]  = "c12";
  varnames [2]  = "c13";
  varnames [3]  = "c14";
  varnames [4]  = "c15";
  varnames [5]  = "c16";
  varnames [6]  = "c22";
  varnames [7]  = "c23";
  varnames [8]  = "c24";
  varnames [9]  = "c25";
  varnames [10] = "c26";
  varnames [11] = "c33";
  varnames [12] = "c34";
  varnames [13] = "c35";
  varnames [14] = "c36";
  varnames [15] = "c44";
  varnames [16] = "c45";
  varnames [17] = "c46";
  varnames [18] = "c55";
  varnames [19] = "c56";
  varnames [20] = "c66";
  varnames [21] = "rho";
      
  ier = ex_get_init ( idexo, cstr, &ndim, &nump, &numel, &numelblk, &numnps, 
    &numess );

  ier = ex_put_init ( idexo, "Title", ndim, nump, numel, numelblk, 
    numnps, numess);
          
  ier = ex_put_var_param ( idexo, "n", 21 );
  
  // TODO Figure out why ier gives (-1) on ex_put_init.  
  
  ier = ex_put_var_names ( idexo, "n", 21, 
    const_cast <char**> ( varnames ));
    
  ier = ex_put_nodal_var ( idexo, 1, 1,  msh.num_nodes, msh.c11 );
  ier = ex_put_nodal_var ( idexo, 1, 2,  msh.num_nodes, msh.c12 );
  ier = ex_put_nodal_var ( idexo, 1, 3,  msh.num_nodes, msh.c13 );
  ier = ex_put_nodal_var ( idexo, 1, 4,  msh.num_nodes, msh.c14 );
  ier = ex_put_nodal_var ( idexo, 1, 5,  msh.num_nodes, msh.c15 );
  ier = ex_put_nodal_var ( idexo, 1, 6,  msh.num_nodes, msh.c16 );
  ier = ex_put_nodal_var ( idexo, 1, 7,  msh.num_nodes, msh.c22 );
  ier = ex_put_nodal_var ( idexo, 1, 8,  msh.num_nodes, msh.c23 );
  ier = ex_put_nodal_var ( idexo, 1, 9,  msh.num_nodes, msh.c24 );
  ier = ex_put_nodal_var ( idexo, 1, 10, msh.num_nodes, msh.c25 );
  ier = ex_put_nodal_var ( idexo, 1, 11, msh.num_nodes, msh.c26 );
  ier = ex_put_nodal_var ( idexo, 1, 12, msh.num_nodes, msh.c33 );
  ier = ex_put_nodal_var ( idexo, 1, 13, msh.num_nodes, msh.c34 );
  ier = ex_put_nodal_var ( idexo, 1, 14, msh.num_nodes, msh.c35 );
  ier = ex_put_nodal_var ( idexo, 1, 15, msh.num_nodes, msh.c36 );
  ier = ex_put_nodal_var ( idexo, 1, 16, msh.num_nodes, msh.c44 );
  ier = ex_put_nodal_var ( idexo, 1, 17, msh.num_nodes, msh.c45 );
  ier = ex_put_nodal_var ( idexo, 1, 18, msh.num_nodes, msh.c46 );
  ier = ex_put_nodal_var ( idexo, 1, 19, msh.num_nodes, msh.c55 );
  ier = ex_put_nodal_var ( idexo, 1, 20, msh.num_nodes, msh.c56 );
  ier = ex_put_nodal_var ( idexo, 1, 21, msh.num_nodes, msh.c66 );
  ier = ex_put_nodal_var ( idexo, 1, 22, msh.num_nodes, msh.rho );
    
}