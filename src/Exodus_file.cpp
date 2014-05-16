#include <stdio.h>
#include <cstdio>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include "classes.hpp"

using namespace std;

void Exodus_file::writeNew ( Mesh &msh )
{
  
  cout << "creating file." << endl;
    
  int id = ex_create ( "./text.ex2", EX_CLOBBER, &comp_ws, &io_ws );
  
  ier = ex_put_init ( id, "TEST", 3, msh.num_nodes, msh.num_elem, 
    msh.refineElemConn.size(), 0, 0 );
  
  ier = ex_put_coord        ( id, msh.xmsh, msh.ymsh, msh.zmsh );
  
  ier = ex_put_node_num_map ( id, msh.node_num_map );
  ier = ex_put_elem_num_map ( id, msh.elem_num_map );
  ier = ex_put_map          ( id, msh.elem_num_map );
  
  int blkNum = 0;
  for ( vector < vector <int> > :: iterator out=msh.refineElemConn.begin(); 
    out!=msh.refineElemConn.end(); ++out )
  {
    
    blkNum++;
    int numElemInBlk = out->size() / double (msh.num_node_per_elem);
    int *blkElemConn = new int 
      [out->size()*msh.num_node_per_elem];

    int connInd = 0;
    for ( vector <int> :: iterator in=out->begin(); in!=out->end(); ++in )
    {
      
      blkElemConn[connInd] = *in;
      connInd++;
      
    }        
    
    ier = ex_put_elem_block   ( id, blkNum, "TETRA", numElemInBlk, 4, 0 );
    ier = ex_put_elem_conn    ( id, blkNum, blkElemConn );
    delete [] blkElemConn;
    
  }

  ier = ex_close ( id );
  
  if ( ier == 0 )
    cout << "New exodus file closed succesfully" << endl;
  
}

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

void Exodus_file::merge ( Region &reg, Model_file &mod ) 
{
  
  cout << "Merging model.\n";
  
  if ( allFiles == true )
  {
    mod.colReg1 = true;
    mod.colReg2 = true;
    mod.lonReg1 = true;
    mod.lonReg2 = true;
    mod.lonReg3 = true;
    mod.lonReg4 = true;
    
    if ( mod.intentions == "CRUST" || mod.intentions == "TOPOGRAPHY" )
    {
      mod.radMin = 6272;
      mod.radReg9 = true;
      mod.radReg8 = true;
      mod.radReg7 = true;
    } 
    // else
    // {
    //   mod.radMin = 0;
    // }
    
  }
      
  if ( mod.colReg1 == true )
    colReg.push_back ("col000-090");
  if ( mod.colReg2 == true )
    colReg.push_back ("col090-180");
    
  if ( mod.lonReg1 == true )
    lonReg.push_back ("lon000-090");
  if ( mod.lonReg2 == true )
    lonReg.push_back ("lon090-180");
  if ( mod.lonReg3 == true )
    lonReg.push_back ("lon180-270");
  if ( mod.lonReg4 == true )
    lonReg.push_back ("lon270-360");
                    
  if ( mod.radReg1 == true )
    radReg.push_back ( "rad0000-1221" );    
  if ( mod.radReg2 == true )
    radReg.push_back ( "rad1221-3480" );  
  if ( mod.radReg3 == true )
    radReg.push_back ( "rad3480-5371" );  
  if ( mod.radReg4 == true )
    radReg.push_back ( "rad5371-5701" );
  if ( mod.radReg5 == true )
    radReg.push_back ( "rad5701-5971" );  
  if ( mod.radReg6 == true )
    radReg.push_back ( "rad5971-6271" );           
  if ( mod.radReg7 == true )
    radReg.push_back ( "rad6271-6319" );            
  if ( mod.radReg8 == true )
    radReg.push_back ( "rad6319-6351" );  
  if ( mod.radReg9 == true )
    radReg.push_back ( "rad6351-6371" );
    
  int l = 0;
  Exodus_file currentExo;
  for ( vector <string>::iterator i=colReg.begin(); i!=colReg.end(); ++i ) 
  {
    for ( vector <string>::iterator j=lonReg.begin(); j!=lonReg.end(); ++j ) 
    {
      
      string dir  = mod.mesh_directory;
      string call = mod.mesh_directory;
      
      dir.append (*i);
      dir.append (".");
      dir.append (*j);
      dir.append ("/");
      
      vector <string> fnames;
      
      for ( vector <string>::iterator l=radReg.begin(); l!=radReg.end(); ++l )
      {
        fnames.push_back ( dir + *i + "." + *j + "." + *l + ".000.ex2" );
      }      
                  
      for ( vector <string>::iterator k=fnames.begin(); k!=fnames.end(); ++k ) 
      {                                
                        
        reg.regionsExo.push_back ( currentExo );
        reg.regionsExo[l].fname = *k;
        l++;
                        
      }
    }
  }  
             
}

void Exodus_file::writeParams ( Mesh &msh )
{

  char *cstr = new char [MAX_LINE_LENGTH];
  const char *varnames[28];

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
  varnames [22] = "Q__";
  varnames [23] = "elv";
  varnames [24] = "siz";
  varnames [25] = "du1";
  varnames [26] = "du2";
  varnames [27] = "du3";
  
  ier = ex_get_init ( idexo, cstr, &ndim, &nump, &numel, &numelblk, &numnps, 
    &numess );

  ier = ex_put_init ( idexo, "Title", ndim, nump, numel, numelblk, 
    numnps, numess);
          
  ier = ex_put_var_param ( idexo, "n", 28 );
  
  // TODO Figure out why ier gives (-1) on ex_put_init.  
  
  ier = ex_put_var_names ( idexo, "n", 28, 
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
  ier = ex_put_nodal_var ( idexo, 1, 23, msh.num_nodes, msh.Q__ );
  ier = ex_put_nodal_var ( idexo, 1, 24, msh.num_nodes, msh.elv );
  ier = ex_put_nodal_var ( idexo, 1, 25, msh.num_nodes, msh.siz );
  ier = ex_put_nodal_var ( idexo, 1, 26, msh.num_nodes, msh.du1 );
  ier = ex_put_nodal_var ( idexo, 1, 27, msh.num_nodes, msh.du2 );
  ier = ex_put_nodal_var ( idexo, 1, 28, msh.num_nodes, msh.du3 );
    
}
