#include "classes.h"
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
    exo.openFile   ( drv.params[i] );
    mshLoc.getInfo ( exo.idexo );
    exo.closeFile  ();    
    
    mshGlb.num_nodes = mshLoc.num_nodes + mshGlb.num_nodes;
    mshGlb.num_elem  = mshLoc.num_elem  + mshGlb.num_elem;
  }
  
  // ********************************************************************* //
  //                            READ IN PARAMS                             //
  // ********************************************************************* //
    
  // Read other model parameters.
  // cout << drv.params[exo.num_mesh_files+5] << std::flush;
  mod.input_model_directory = drv.params[exo.num_mesh_files+1];
  mod.input_model_file_type = drv.params[exo.num_mesh_files+2];
  mod.input_model_physics   = drv.params[exo.num_mesh_files+3];
  mod.absolute_or_perturb   = drv.params[exo.num_mesh_files+4]; 
  mod.rotAng                = stod (drv.params[exo.num_mesh_files+5]);
  mod.rotVecX               = stod (drv.params[exo.num_mesh_files+6]);
  mod.rotVecY               = stod (drv.params[exo.num_mesh_files+7]);
  mod.rotVecZ               = stod (drv.params[exo.num_mesh_files+8]);
  mod.intentions            = drv.params[exo.num_mesh_files+9];
  mod.dimensions            = drv.params[exo.num_mesh_files+10];
  mod.interpolation         = drv.params[exo.num_mesh_files+11];
       
  cout << "\nModel information:\n* INPUT_MODEL_DIRECTORY: "   <<
    mod.input_model_directory << "\n* INPUT_MODEL_FILE_TYPE: " <<
    mod.input_model_file_type << "\n* INPUT_MODEL_PHYSICS: "   <<
    mod.input_model_physics   << "\n* ABSOLUTE_OR_PERTURB: "   <<
    mod.absolute_or_perturb   << "\n";
  
  // ********************************************************************* //
  //                            READ IN MODEL                              //
  // ********************************************************************* //
  
  mod.read ();  
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
    
    exo.openFile         ( drv.params[i] );
    mshLoc.getInfo       ( exo.idexo );
    mshLoc.allocateMesh  ( mshLoc.num_nodes );
    mshLoc.populateCoord ( exo.idexo ); 
    
    cout << "Changing physics.\n";
    mod.openUp ( mshLoc );
    cout << "Physics has been changed.\n";
    
    if ( mod.intentions != "INITIALIZE" ) {
      mshLoc.populateParams ( exo.idexo );
    }
    
    if ( mod.dimensions == "NORMALIZED" ) {
      for ( i=0; i<mshLoc.num_nodes; i++ ) {
        mshLoc.xmsh[i] = mshLoc.xmsh[i] * con.R_EARTH;
        mshLoc.ymsh[i] = mshLoc.ymsh[i] * con.R_EARTH;
        mshLoc.zmsh[i] = mshLoc.zmsh[i] * con.R_EARTH;
      }
    }
         
    mshLoc.interpolateModel ( mod );
    exo.closeFile           ();        
    
  }  
  return 0;
  
}