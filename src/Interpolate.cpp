#include "classes.hpp"

// using namespace std;

int main () 
{
  
  std::ifstream myfile;
  
  Region        reg;  
  Constants     con;
  Exodus_file   exo;
  Driver        drv;
  Model_file    mod;
  Utilities     utl;
  Discontinuity dis;
  
  std::cout << "Begin model building.\n";  
      
  // Read parameter file.
  drv.initialize ( mod, dis, utl, exo, reg );  
  
  if ( mod.intentions != "INTERPOLATE" ) 
  {
    
    std::cout << "\nHey bru I think you've got a mistake in your parameter " << 
      "file. I'm stopping here for you to check it out ( wrong intention? )\n" 
        << std::endl;
    exit ( EXIT_FAILURE );
  }
  
  
  if ( mod.intentions == "INTERPOLATE" )   
  {
    
    std::cout << "\n----- Interpolating -----\n";
    std::cout << "\n";
    
    mod.createKDTreeUnpacked ( );
    dis.createKDTreePacked   ( );
    
    for ( std::vector < Exodus_file > :: 
      iterator exoFile=reg.regionsExo.begin();
      exoFile!=reg.regionsExo.end(); ++exoFile )       
    {  
        
      std::cout << "\n";

      Mesh msh;
      Interpolator ipl;    
    
      exoFile -> openFile    ( exoFile -> fname );        
      msh.getInfo            ( exoFile -> idexo, 'p' );         
      ipl.interpolate        ( msh, mod, dis );
      exoFile -> writeParams ( msh );
      exoFile -> closeFile   ( );      
    
    }
  }
                  
  return 0;    
}