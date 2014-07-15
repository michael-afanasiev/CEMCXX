#include "classes.hpp"

int main () 
{
  
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
  drv.checkUsage ( mod, "TOPOGRAPHY" );
  
  // Read topography and whatnot.   
  dis.readTopography       ( );  
  dis.createKDTreeUnpacked ( );
      
  std::cout << "\n----- Interpolating Topography -----\n";
  std::cout << "\n";
    
  for ( std::vector < Exodus_file > :: 
    iterator exoFile=reg.regionsExo.begin();
    exoFile!=reg.regionsExo.end(); ++exoFile )       
  {  
      
    std::cout << "\n";

    Mesh msh;
    Interpolator ipl;    
  
    exoFile -> openFile    ( exoFile -> fname );        
    msh.getInfo            ( exoFile -> idexo );    
    ipl.interpolateTopo    ( msh, dis );
    exoFile -> writeParams ( msh );
    msh.deallocateMesh     ( mod );
    exoFile -> closeFile   ( );      
  
  }
                  
  return 0;    
}
