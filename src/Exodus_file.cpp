#include <stdio.h>
#include <cstdio>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include "classes.hpp"

using namespace std;

void Exodus_file::openFile ( string fname ) 
{
  /**
   * Opens an exodus file and populates the idexo field in the exodus 
   * object.
   */
    
  std::cout << "Opening exodus file: " << fname << std::flush << std::endl;
    
  idexo = ex_open ( fname.c_str(), EX_WRITE, &comp_ws, &io_ws, &vers ); 
  
  if (idexo < 0) 
  {
    std::cout << "***Fatal error opening exodus file. Exiting." << std::endl;  
    exit (EXIT_FAILURE);  
  } 
  
} 

void Exodus_file::closeFile ( ) 
{
  /**
   * Closes the exodus file associated with the object.
   */
  
  ier = ex_close ( idexo );
  
  if (ier == 0) 
  {
    std::cout << "File closed succesfully." << std::flush << std::endl;
  } 
  else 
  {
    std::cout << "***Fatal error closing exodus file. Exiting" << std::endl;
    exit (EXIT_FAILURE);  
  }
}

void Exodus_file::merge ( Region &reg, Model_file &mod ) 
{
  /**
   * This function helps decide which parts of the cem to open, given
   * a model request. The actual mod.---Reg fields are populated in the
   * utilites. rotate function. That could probably be made cleaner. We 
   * basically get 'regionsExo' out of this, which is the master vector of
   * all exodus files are needed in this particular run.
   */
    
  std::cout << "Merging model." << std::flush << std::endl;
  
  /* An option to include the entire globe. This is useful for things that
   * don't necessary have a minimum radius or lat/lon bounds (like the crust,
   * or topography). I chose 100 km as where to stop the crust and topography
   * interpolation, because there's no deeper crust than 100 km. Allfiles is
   * automatically set if you choose crust or topography. */ 
  if ( allFiles == true )
  {
    mod.colReg1 = true;
    mod.colReg2 = true;
    mod.lonReg1 = true;
    mod.lonReg2 = true;
    mod.lonReg3 = true;
    mod.lonReg4 = true;
    
    // These files define the top 100 km.
    if ( mod.intentions == "CRUST" || mod.intentions == "TOPOGRAPHY" )
    {
      mod.radMin   = 6271.;
      mod.radReg32 = true;
      mod.radReg33 = true;
      mod.radReg34 = true;
      mod.radReg35 = true;
      mod.radReg36 = true;
      mod.radReg37 = true;
    } 
        
  }
   
  /* Decide which regions to add to a normal interpolation / extraction. These
   * values are set in utilities.rotate, and are pulled from the maximum and 
   * minimum xyz values read in from the model. */   
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
       
  if ( mod.radReg1   == true )
    radReg.push_back ( "rad0000-1221" );        
  if ( mod.radReg2   == true )
    radReg.push_back ( "rad1221-3480" );                
  if ( mod.radReg3   == true )
    radReg.push_back ( "rad3480-3571" );              
  if ( mod.radReg4   == true )
    radReg.push_back ( "rad3571-3671" );              
  if ( mod.radReg5   == true )
    radReg.push_back ( "rad3671-3771" );              
  if ( mod.radReg6   == true )
    radReg.push_back ( "rad3771-3871" );              
  if ( mod.radReg7   == true )
    radReg.push_back ( "rad3871-3971" );              
  if ( mod.radReg8   == true )
    radReg.push_back ( "rad3971-4071" );              
  if ( mod.radReg9   == true )
    radReg.push_back ( "rad4071-4171" );              
  if ( mod.radReg10   == true )
    radReg.push_back ( "rad4171-4271" );              
  if ( mod.radReg11   == true )
    radReg.push_back ( "rad4271-4371" );              
  if ( mod.radReg12   == true )
    radReg.push_back ( "rad4371-4471" );              
  if ( mod.radReg13   == true )
    radReg.push_back ( "rad4471-4571" );              
  if ( mod.radReg14   == true )
    radReg.push_back ( "rad4571-4671" );              
  if ( mod.radReg15   == true )
    radReg.push_back ( "rad4671-4771" );              
  if ( mod.radReg16   == true )
    radReg.push_back ( "rad4771-4871" );              
  if ( mod.radReg17   == true )
    radReg.push_back ( "rad4871-4971" );              
  if ( mod.radReg18  == true )
    radReg.push_back ( "rad4971-5071" );              
  if ( mod.radReg19  == true )
    radReg.push_back ( "rad5071-5171" );              
  if ( mod.radReg20  == true )
    radReg.push_back ( "rad5171-5271" );              
  if ( mod.radReg21  == true )
    radReg.push_back ( "rad5271-5371" );              
  if ( mod.radReg22  == true )
    radReg.push_back ( "rad5371-5426" );            
  if ( mod.radReg23  == true )
    radReg.push_back ( "rad5426-5481" );              
  if ( mod.radReg24  == true )
    radReg.push_back ( "rad5481-5536" );              
  if ( mod.radReg25  == true )
    radReg.push_back ( "rad5536-5591" );              
  if ( mod.radReg26  == true )
    radReg.push_back ( "rad5591-5646" );              
  if ( mod.radReg27  == true )
    radReg.push_back ( "rad5646-5701" );              
  if ( mod.radReg28 == true )
    radReg.push_back ( "rad5701-5746" );              
  if ( mod.radReg29 == true )
    radReg.push_back ( "rad5746-5791" );              
  if ( mod.radReg30 == true )
    radReg.push_back ( "rad5791-5836" );           
  if ( mod.radReg31 == true )
    radReg.push_back ( "rad5836-5881" );             
  if ( mod.radReg32 == true )
    radReg.push_back ( "rad5881-5926" );              
  if ( mod.radReg33 == true )
    radReg.push_back ( "rad5926-5971" );             
  if ( mod.radReg34 == true )
    radReg.push_back ( "rad5971-6021" );              
  if ( mod.radReg35 == true )
    radReg.push_back ( "rad6021-6071" );              
  if ( mod.radReg36 == true )
    radReg.push_back ( "rad6071-6121" );              
  if ( mod.radReg37 == true )
    radReg.push_back ( "rad6121-6171" );              
  if ( mod.radReg38 == true )
    radReg.push_back ( "rad6171-6221" );              
  if ( mod.radReg39 == true )
    radReg.push_back ( "rad6221-6271" );              
  if ( mod.radReg40 == true )
    radReg.push_back ( "rad6271-6291" );            
  if ( mod.radReg41 == true )
    radReg.push_back ( "rad6291-6311" );              
  if ( mod.radReg42 == true )
    radReg.push_back ( "rad6311-6331" );              
  if ( mod.radReg43 == true )
    radReg.push_back ( "rad6331-6351" );              
  if ( mod.radReg44 == true )
    radReg.push_back ( "rad6351-6361" );              
  if ( mod.radReg45 == true )
    radReg.push_back ( "rad6361-6371" );
    
  int l      = 0;  
  string dir = mod.mesh_directory;
  Exodus_file currentExo;
  /* Loop over col/lon regions. Within each region we will build radially. */
  for ( vector <string>::iterator i=colReg.begin(); i!=colReg.end(); ++i ) 
  {
    for ( vector <string>::iterator j=lonReg.begin(); j!=lonReg.end(); ++j ) 
    {
                  
      /* We build up a vector of appropriate radial regions for each col/lon
       * chunk. */
      vector <string> fnames;
      for ( vector <string>::iterator ll=radReg.begin(); ll!=radReg.end(); 
        ++ll )
      {
        fnames.push_back ( dir + *i + "." + *j + "." + *ll + ".000.ex2" );
        reg.colReg.push_back ( *i );
        reg.lonReg.push_back ( *j );
      }      
                  
      /* Push back the details of the exodus file in currentExo. As we do this,
       * put the name of the file in the region object as well. In the end we
       * have the file details and filename stored in the Exodus_file object. */
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
  
 /**
   * This guy writes parameters back to the exodus files, for an interpolation
   * step.
   */
  

  char *cstr = new char [MAX_LINE_LENGTH];
  const char *varnames[28];

  int ndim;
  int nump;
  int numel;
  int numelblk;
  int numnps;
  int numess;
  
  // Parameter names.
  varnames[0]  = "c11";
  varnames[1]  = "c12";
  varnames[2]  = "c13";
  varnames[3]  = "c14";
  varnames[4]  = "c15";
  varnames[5]  = "c16";
  varnames[6]  = "c22";
  varnames[7]  = "c23";
  varnames[8]  = "c24";
  varnames[9]  = "c25";
  varnames[10] = "c26";
  varnames[11] = "c33";
  varnames[12] = "c34";
  varnames[13] = "c35";
  varnames[14] = "c36";
  varnames[15] = "c44";
  varnames[16] = "c45";
  varnames[17] = "c46";
  varnames[18] = "c55";
  varnames[19] = "c56";
  varnames[20] = "c66";
  varnames[21] = "rho";
  varnames[22] = "Q__";
  varnames[23] = "elv";
  varnames[24] = "siz";
  varnames[25] = "du1";
  varnames[26] = "du2";
  varnames[27] = "du3";
  
  /* Get the current parameters from the open file. We will basically re-write
   * the same parameters back to the original file, with updated tensor 
   * components */
  ier = ex_get_init ( idexo, cstr, &ndim, &nump, &numel, &numelblk, &numnps, 
    &numess );

//  ier = ex_put_init ( idexo, "Title", ndim, nump, numel, numelblk, 
//    numnps, numess);
          
  ier = ex_put_var_param ( idexo, "n", 28 );
  
//  ier = ex_put_var_names ( idexo, "n", 28, 
//    const_cast <char**> ( varnames ));
    
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
