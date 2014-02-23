#include <stdio.h>
#include <cstdio>
#include <string.h>
#include <unistd.h>
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

void Exodus_file::merge ( Region &reg, Model_file &mod ) 
{
  
  cout << "Merging model.\n";
      
  if ( mod.colMax <= 90. ) {
    colReg.push_back ("col0-90");    
  }
  if ( (mod.colMax >= 90.) && (mod.colMin >= 90.) ) {
    colReg.push_back ("col90-180");
  }
  if ( (mod.colMax >= 90.) && (mod.colMin <= 90.) ) {
    colReg.push_back ("col0-90");
    colReg.push_back ("col90-180");
  } 
  
  double lonMinLoc = mod.lonMin;
  double lonMaxLoc = mod.lonMax;
  
  for ( int l=0; l<=270; l+=90 ) {
      
    if ( lonMinLoc < 0. ) {
      lonMinLoc = 180. - lonMinLoc;                
    }
    if ( lonMaxLoc < 0. ) {
      lonMaxLoc = 180. - lonMaxLoc;        
    }
    
    if ( (lonMinLoc >= l) && (lonMaxLoc <= l+90) ) {
      
      string dum1 = to_string (l);
      string dum2 = to_string (l+90);            
      string full = "lon";
      
      full.append (dum1);
      full.append ("-");
      full.append (dum2);
      
      lonReg.push_back (full);      
            
    }      
  }
    
  if ( mod.radMin <= 1221 ) {
    radReg.push_back ( "rad0-1221" );
  }
    
  if ( mod.radMin <= 3480 ) {
    radReg.push_back ( "rad1221-3480" );
  }
    
  if ( mod.radMin <= 5371 ) {
    radReg.push_back ( "rad3480-5371" );
  }
  
  if ( mod.radMin <= 6271 ) {        
    radReg.push_back ( "rad5371-6271" );        
  }
   
  if ( mod.radMin <= 6319 ) {
    radReg.push_back ( "rad6271-6319" );        
  }
    
  if ( mod.radMin <= 6351 ) {
    radReg.push_back ( "rad6319-6351" );
  }
   
  if ( mod.radMin <= 6371 ) {    
    radReg.push_back ( "rad6351-6371" );
  }
    
  int   num;
  float dum1;
  char  dum2;
  int l = 0;
  Exodus_file exii;
  string masterCall = "ejoin -output ./dat/input.ex2 ";
  for ( vector <string>::iterator i=colReg.begin(); i!=colReg.end(); ++i ) {
    for ( vector <string>::iterator j=lonReg.begin(); j!=lonReg.end(); ++j ) {
      for ( vector <string>::iterator k=radReg.begin(); k!=radReg.end(); ++k ) {
        
        string call = mod.mesh_directory;
                        
        call.append (*i);
        call.append (".");
        call.append (*j);
        call.append (".");
        call.append (*k);
        call.append (".ex2");
                
        masterCall.append (call);
        masterCall.append (" ");
        
        reg.regionsExo.push_back ( exii );
        reg.regionsExo[l].fname = call;
        
        // Get number of blocks for later destruction.
        idexo = ex_open    ( call.c_str(), EX_READ, &comp_ws, &io_ws, &vers );
        ier   = ex_inquire ( idexo, EX_INQ_ELEM_BLK, &num, &dum1, 
          &dum2 );        
        ier   = ex_close   ( idexo );
        
        totalBlocks.push_back ( num );
                
        l++;
                
      }
    }
  }  

 
  // system ( masterCall.c_str() );  
             
}

void Exodus_file::splitBack ()
{

  int fBlock = 0;
  int iBlock = 0;
  for ( vector <string>::iterator i=colReg.begin(); i!=colReg.end(); ++i ) {
    for ( vector <string>::iterator j=lonReg.begin(); j!=lonReg.end(); ++j ) {
      for ( vector <string>::iterator k=radReg.begin(); k!=radReg.end(); ++k ) {
        
        
        ofstream myfile ( "./split.txt", ios::out );
        myfile << "delete block all\nundelete block " << 
          fBlock + 1 << " to " << totalBlocks[iBlock] + fBlock << "\nexit";
        
        string call = "cat split.txt | grepos ./dat/input.ex2 ./dat/";
        
        call.append (*i);
        call.append (".");
        call.append (*j);
        call.append (".");
        call.append (*k);
        call.append (".ex2");        
        
        cout << "Running grepos. Splitting region: " << iBlock << "\n" << 
          std::flush;
        myfile.close();
        
        FILE *fp = NULL; 
        fp = popen ( call.c_str(), "r" );
        sleep ( 20 );        
        pclose (fp);
        
        fBlock += totalBlocks[iBlock];
        iBlock++;
                              
      }
    }
  }
}

void Exodus_file::writeParams ( Mesh &msh )
{

  char *cstr = new char [MAX_LINE_LENGTH];
  const char *varnames[27];

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
  varnames [24] = "du1";
  varnames [25] = "du2";
  varnames [26] = "du3";
  
  ier = ex_get_init ( idexo, cstr, &ndim, &nump, &numel, &numelblk, &numnps, 
    &numess );

  ier = ex_put_init ( idexo, "Title", ndim, nump, numel, numelblk, 
    numnps, numess);
          
  ier = ex_put_var_param ( idexo, "n", 27 );
  
  // TODO Figure out why ier gives (-1) on ex_put_init.  
  
  ier = ex_put_var_names ( idexo, "n", 27, 
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
  ier = ex_put_nodal_var ( idexo, 1, 25, msh.num_nodes, msh.du1 );
  ier = ex_put_nodal_var ( idexo, 1, 26, msh.num_nodes, msh.du2 );
  ier = ex_put_nodal_var ( idexo, 1, 27, msh.num_nodes, msh.du3 );
    
}


void Exodus_file::writeSize ( Mesh &msh )
{
  
  char *cstr = new char [MAX_LINE_LENGTH];
  const char *varnames[1];

  int ndim;
  int nump;
  int numel;
  int numelblk;
  int numnps;
  int numess;
  
  varnames [0] = "siz";
  
  ier = ex_get_init ( idexo, cstr, &ndim, &nump, &numel, &numelblk, &numnps, 
    &numess );

  ier = ex_put_init ( idexo, "Title", ndim, nump, numel, numelblk, 
    numnps, numess);
          
  ier = ex_put_var_param ( idexo, "n", 1 );
  
  // TODO Figure out why ier gives (-1) on ex_put_init.  
  
  ier = ex_put_var_names ( idexo, "n", 1, 
    const_cast <char**> ( varnames ));
    
  ier = ex_put_nodal_var ( idexo, 1, 1,  msh.num_nodes, msh.siz );
  
}