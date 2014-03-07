#include "classes.hpp"

// using namespace std;

int main () 
{
  
  std::ofstream myfile;
  
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
  
  if ( mod.intentions != "REFINE" ) 
  {
    
    std::cout << "\nHey bru I think you've got a mistake in your parameter " << 
      "file. I'm stopping here for you to check it out ( wrong intention? )\n" 
        << std::endl;
    exit ( EXIT_FAILURE );
  }
  
  int i = 0;
  if ( mod.intentions == "REFINE" ) 
  {
    
    std::cout << "\n----- Outputting refinment. -----\n";
        
    for ( std::vector < Exodus_file > :: 
      iterator exoFile=reg.regionsExo.begin();
      exoFile!=reg.regionsExo.end(); ++exoFile ) 
    {
        
      std::string fname ( "./dat/refine/refine" );
      fname += std::to_string (i);
      fname += ".txt";
      
      myfile.open ( fname, std::ios::out );        
      std::cout << "\n";        
    
      Mesh         msh;
      Interpolator ipl;
    
      exoFile -> openFile        ( exoFile -> fname );        
      msh.getInfo                ( exoFile -> idexo, 'p' );    
      msh.getNodeNumMap          ( exoFile -> idexo );
      msh.getElemNumMap          ( exoFile -> idexo );
      msh.getElementConnectivity ( exoFile -> idexo );
      ipl.findNodes              ( msh, mod );
      msh.deallocateMesh         ( mod );
      exoFile -> closeFile       ( );  
      
      myfile << exoFile -> fname << std::endl;
      for ( std::vector < int > :: iterator node=ipl.refineArr.begin();
        node!=ipl.refineArr.end(); ++node )
      {
        myfile << *node << std::endl;
      }
      
      myfile.close();
      i++;
    }    
  }
  
  return 0;
}