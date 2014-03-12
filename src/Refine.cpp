#include "classes.hpp"

// using namespace std;

int main ( int argc, char* argv[] ) 
{
  
  std::ofstream myfile;
  
  Region        reg;  
  Constants     con;
  Exodus_file   exo;
  Driver        drv;
  Model_file    mod;
  Utilities     utl;
  Discontinuity dis;
  
  if ( argc < 2 )
  {
    std::cout << "Incorrect usage." << std::endl;
    return 0;
  }
  
  std::string arg = argv[1];
  mod.refineSize = stoi ( arg );
  std::cout << mod.refineSize << std::endl;
  
  std::cout << "Begin model building.\n";
  
  // Read parameter file.
  drv.initialize ( mod, dis, utl, exo, reg );  
  
  if ( mod.intentions != "REFINE" ) 
  {
    
    std::cout << "\nHey bru I think you've got a mistake in your parameter " 
      << "file. I'm stopping here for you to check it out "
      << "( wrong intention? )" << std::endl;
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
      
      std::cout << "\n";        
    
      Mesh         msh;
      Interpolator ipl;
    
      exoFile -> openFile        ( exoFile -> fname );        
      msh.getInfo                ( exoFile -> idexo, 'p' );    
      ipl.findNodes              ( msh, mod );
      exoFile -> writeParams     ( msh );
      msh.deallocateMesh         ( mod );
      exoFile -> closeFile       ( );   
    }    
  }
  
  return 0;
}