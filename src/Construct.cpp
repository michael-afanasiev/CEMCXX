#include "classes.hpp"
#include "kdtree.h"

using namespace std;

int main () 
{
  
  int i;
  
  ifstream myfile;
  
  Constants    con;
  Exodus_file  exo;
  Driver       drv;
  Mesh         msh;
  Model_file   mod;
  Utilities    util;
  Interpolator inter;
  kdtree       *tree;
  
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
    
  // ********************************************************************* //
  //                            QUERY THE MESH                             //
  // ********************************************************************* //

  exo.openFile ( drv.params[0] );
  msh.getInfo  ( exo.idexo );
  
  // ********************************************************************* //
  //                            READ IN PARAMS                             //
  // ********************************************************************* //
  
  mod.populateParams ( drv, exo );

  // ********************************************************************* //
  //                            READ IN MODEL                              //
  // ********************************************************************* //
  
  mod.read ();  
  
  // if ( mod.rotAng != 0 ) {
  //   mod.rotRad = mod.rotAng * con.PI / con.o80;   
  //   util.rotate ( mod );
  // }    
        
  // ********************************************************************* //
  //                            INTERPOLATE                                //
  // ********************************************************************* //
      
  // Report to the friendly user.    
  cout << "\nWe're working with " << msh.num_nodes << " nodes, and " << 
    msh.num_elem << " elements.\n";
                          
  msh.getInfo        ( exo.idexo );
  msh.allocateMesh   ( msh.num_nodes );    
  
  msh.populateCoord  ( exo.idexo ); 
  msh.populateParams ( exo.idexo, mod );
                       
  mod.openUp         ();  // Dome rock.. dome rock.. dome rock
  
  msh.reNormalize    ( mod );    
  
  // inter.interpolate  ( msh, mod );
  inter.exterpolator ( msh, exo, mod );
  
  exo.writeParams    ( msh );
  exo.closeFile      ();        
                   
  return 0;
  
}