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
  
  std::cout << "Begin model building." << std::endl;  
      
  // Read parameter file.
  drv.initialize ( mod, dis, utl, exo, reg );  
  drv.checkUsage ( mod, "INTERPOLATE" );
      
  std::cout << "\n----- Interpolating -----\n";
  std::cout << "\n";
  
  dis.createKDTreePacked ( );
  
  for ( std::vector < Exodus_file > :: 
    iterator exoFile=reg.regionsExo.begin();
    exoFile!=reg.regionsExo.end(); ++exoFile )       
  {  
      
    std::cout << "\n";

    Mesh msh;
    Interpolator ipl;    
  
    exoFile -> openFile      ( exoFile -> fname );        
    msh.getInfo              ( exoFile -> idexo, 'p' );  
    mod.createKDTreeUnpacked ( msh );       
    ipl.interpolate          ( msh, mod, dis );
    exoFile -> writeParams   ( msh );
    exoFile -> closeFile     ( );     
    msh.deallocateMesh       ( mod );
    mod.deallocate           ( );
    
  }
  
  drv.report ( mod );              
  return 0;    
}
