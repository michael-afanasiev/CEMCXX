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
      
  exo.allFiles = true;
      
  // Read parameter file.
  drv.initialize ( mod, dis, utl, exo, reg );  
    
  if ( mod.intentions != "CRUST" ) 
  {
    
    std::cout << "\nHey bru I think you've got a mistake in your parameter " << 
      "file. I'm stopping here for you to check it out ( wrong intention? )\n" 
        << std::endl;
    exit ( EXIT_FAILURE );
  }
  
  
  if ( mod.intentions == "CRUST" )   
  {
    
    std::cout << "\n----- Interpolating Crust -----\n";
    std::cout << "\n";
    
    dis.createKDTreePacked   ( );
    
    for ( std::vector < Exodus_file > :: 
      iterator exoFile=reg.regionsExo.begin();
      exoFile!=reg.regionsExo.end(); ++exoFile )       
    {  
        
      std::cout << "\n";

      Mesh msh;
      Interpolator ipl;    
    
      exoFile -> openFile    ( exoFile -> fname );        
      msh.getInfo            ( exoFile -> idexo );    
      ipl.interpolateCrust   ( msh, dis, mod );
      exoFile -> writeParams ( msh );
      msh.deallocateMesh     ( mod );
      exoFile -> closeFile   ( );      
    
    }
  }
                  
  return 0;    
}
