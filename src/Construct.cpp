#include "classes.hpp"
#include <time.h>
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
  //                            READ IN PARAMS                             //
  // ********************************************************************* //
  
  mod.populateParams ( drv, exo );
        
  // ********************************************************************* //
  //                            QUERY THE MESH                             //
  // ********************************************************************* //

  // exo.openFile ( drv.params[0] );

  // ********************************************************************* //
  //                            READ IN MODEL                              //
  // ********************************************************************* //

  
  mod.read ();  
  mod.readDiscontinuities ();
  
  // if ( mod.rotAng != 0 ) {
  //   mod.rotRad = mod.rotAng * con.PI / con.o80;   
  //   util.rotate ( mod );
  // }    
  
  // ********************************************************************* //
  //                           CONSTRUCT THE MESH                          //
  // ********************************************************************* //
  
  exo.merge ( mod );
  exo.openFile ( "./dat/input.ex2" );
        
  // ********************************************************************* //
  //                            INTERPOLATE                                //
  // ********************************************************************* //
      
  if ( mod.intentions == "INTERPOLATE" )  {                             
    msh.getInfo        ( exo.idexo, 'p' );
                       
    mod.openUp         ();  // Dome rock.. dome rock.. dome rock
  
    Interpolator inter;
    inter.interpolate  ( msh, mod );
  
    exo.writeParams    ( msh );
    exo.closeFile      ();        
    
    exo.splitBack ();
  }
  // ********************************************************************* //
  //                               EXTRACT                                 //
  // ********************************************************************* //
  
  if ( mod.intentions == "EXTRACT" ) {        
    msh.getInfo        ( exo.idexo, 'p' );
    
    mod.openUp          ();
    
    Interpolator inter;
    inter.exterpolator ( msh, exo, mod );
    
    mod.writeSES3D     ();
  
    exo.closeFile      ();                
  } 
  return 0;
  
}