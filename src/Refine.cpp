#include "classes.hpp"
#include <cmath>
#include <cstring>

void checkArgv ( )
{
  std::cerr << "Usage: refine [desired refinement size]" << std::endl;
  exit ( EXIT_FAILURE );
}

int main ( int argc, char* argv[] ) 
{
 
  std::ofstream      splitFile;
  
  Region        reg;  
  Exodus_file   exo;
  Driver        drv;
  Utilities     utl;
  Model_file    mod;
  Discontinuity dis;
  
  double toMB = 9.5367e-7;
  splitFile.open ( "./splitIntructions.txt", std::ios::out );
  
  if ( argc < 2 )
    checkArgv ( );
  
  std::string arg = argv[1];
  mod.refineSize = stoi ( arg );
  std::cout << mod.refineSize << std::endl;
  
  std::cout << "Begin model building.\n";
  
  // Read parameter file.
  drv.initialize ( mod, dis, utl, exo, reg ); 
  drv.checkUsage ( mod, "REFINE" );
  
  std::cout << "\n----- Outputting refinment. -----\n";
      
  int i = 0;
  for ( std::vector < Exodus_file > :: 
    iterator exoFile=reg.regionsExo.begin();
    exoFile!=reg.regionsExo.end(); ++exoFile ) 
  {            
    std::cout << "\n";        
  
    FILE          *pFile;
    Mesh           msh;
    Interpolator   ipl;
    std::ofstream  myfile;    
            
    std::string fname = "./dat/refine/Refine";
    fname.append ( std::to_string(i) );
    fname.append ( ".txt" );
    
    myfile.open ( fname, std::ios::out );
    
    myfile << exoFile -> fname << std::endl;
  
    exoFile -> openFile ( exoFile -> fname );        
    msh.getInfo         ( exoFile -> idexo, 'p' );
    
    // Get current file size.
    pFile = fopen ( exoFile -> fname.c_str(), "rb" );
    fseek ( pFile, 0, SEEK_END );
    int size = ftell ( pFile );
    fclose ( pFile );        
    
    msh.getNodeNumMap     ( exoFile -> idexo );    
    ipl.findNodes         ( msh, mod, myfile );
        
    double percent = double (msh.numFound) / double (msh.num_nodes);
    std::cout << "Percentage of refined elements: " << percent << std::endl;
    double newElem = msh.num_nodes - msh.numFound + msh.numFound * pow ( 8, 2 );
    double newSize = ( newElem / msh.num_elem ) * size * toMB;
    
    std::cout << "Old size: "           << size * toMB << std::endl;
    std::cout << "Estimated new size: " << newSize << std::endl;
    
    if ( newSize > 250. )
    {
      std::cout << "I think " << exoFile -> fname << " will be larger than "
      << "250 MB after refinment. Splitting file ... " << std::endl;
      splitFile << exoFile -> fname << std::endl;
    }
    
    msh.deallocateMesh     ( mod );
    exoFile -> closeFile   ( );   
    
    i++;
  }    
  
  splitFile.close ();
  return 0;
}