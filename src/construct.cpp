#include "classes.hpp"

using namespace std;

int main () 
{
  
  int i;
  
  ifstream myfile;
  
  Constants   con;
  Exodus_file exo;
  Driver      drv;
  Mesh        mshLoc;
  Mesh        mshGlb;
  Model_file  mod;
  Utilities   util;
  
  cout << "Begin model building.\n";  
  
  // ********************************************************************* //
  //                       READING THE PARAMETER FILE                      //
  // ********************************************************************* //
    
  // Determine how many lines to read from parameter file.
  drv.openDriver ( myfile );
  
  // Read those lines from parameter file and put them into 'params'.
  drv.readDriver ( myfile );
  
  // Close driver.
  drv.closeDriver ( myfile );
  
  // Determine the number of mesh files.
  exo.num_mesh_files = stoi ( drv.params[0] );
  cout << "Reading " << exo.num_mesh_files << " mesh file(s).\n";
    
  // ********************************************************************* //
  //                            QUERY THE MESH                             //
  // ********************************************************************* //
    
  // Since we've got a variable number of mesh files, 
  // the next step is mandatory. We're counting the total number of elements
  // and nodes across the files. This might be able to be cleaned up a teeny
  // bit by initializing seperate Exodus_file object for each open file, but
  // I don't think it's worth it right now.
  mshGlb.num_nodes = 0;
  mshGlb.num_elem  = 0;
  for ( i=1; i<exo.num_mesh_files+1; i++ ) {
    exo.openFile ( drv.params[i] );
    mshLoc.getInfo ( exo.idexo );
    
    mshGlb.num_nodes = mshLoc.num_nodes + mshGlb.num_nodes;
    mshGlb.num_elem  = mshLoc.num_elem  + mshGlb.num_elem;
  }
  
  // ********************************************************************* //
  //                            READ IN PARAMS                             //
  // ********************************************************************* //
  
  mod.populateParams ( drv, exo );

  // ********************************************************************* //
  //                            READ IN MODEL                              //
  // ********************************************************************* //
  
  mod.read ();  
  mod.findMinMax ();
  
  if ( mod.input_model_file_type == "SES3D" ) {    
    mod.populateRadiansSES3D ();
    mod.colLonRad2xyzSES3D   ();    
  }
  
  if ( mod.rotAng != 0 ) {
    mod.rotRad = mod.rotAng * con.PI / con.o80;   
    util.rotate ( mod );
  }    
        
  // ********************************************************************* //
  //                            INTERPOLATE                                //
  // ********************************************************************* //
      
  // Report to the friendly user.    
  cout << "\nWe're working with " << mshGlb.num_nodes << " nodes, and " << 
    mshGlb.num_elem << " elements.\n";
    
  for ( i=1; i<exo.num_mesh_files+1; i++ ) { 
     
    exo.openFile            ( drv.params[i] );
                            
    mshLoc.getInfo          ( exo.idexo );
    mshLoc.allocateMesh     ( mshLoc.num_nodes );
    mshLoc.populateCoord    ( exo.idexo ); 
    mshLoc.populateParams   ( exo.idexo, mod );
                            
    mod.openUp              ();  // Dome rock.. dome rock.. dome rock
    
    mshLoc.reNormalize      ( mod );    
    mshLoc.interpolateModel ( mod );
    
    exo.closeFile           ();        
    
    
    }
               
  return 0;
  
}