#include "classes.hpp"
#include <cmath>
#include <cstring>

void estimateSize ( Mesh &, int &, bool &t, std::ofstream &, std::string  );
void setupRefineFile ( std::ofstream &, std::string, int );

int main ( ) 
{
 
  std::ofstream splitFile;
  
  Region        reg;  
  Exodus_file   exo;
  Driver        drv;
  Utilities     utl;
  Model_file    mod;
  Discontinuity dis;
  Constants     con;
  
  splitFile.open ( "./tmp/splitInstructions.txt", std::ios::out );
  
  std::cout << "Begin model building.\n";
  
  // Read parameter file.
  drv.initialize ( mod, dis, utl, exo, reg );   
  drv.checkUsage ( mod, "REFINE" );
  
  std::cout << "\n----- Outputting refinment. -----\n";
      
  int numRefine  = 0;
  bool willSplit = false;
  
  for ( std::vector < Exodus_file > :: 
    iterator exoFile=reg.regionsExo.begin();
    exoFile!=reg.regionsExo.end(); ++exoFile ) 
  {         
       
    std::cout << "\n";        
  
    Mesh           msh;
    Interpolator   ipl;            
  
    // Initialize the exodus file.
    exoFile -> openFile ( exoFile -> fname );        
    msh.getInfo         ( exoFile -> idexo, 'p' );
    
    // Get current file size.
    int size = utl.getFilesize ( exoFile -> fname );
    
    // Set up file to hold refined elements.
    std::ofstream  myfile;    
    setupRefineFile ( myfile, exoFile -> fname, numRefine );
    
    // Get node number map and find refinable nodes.

    msh.getNodeNumMap          ( exoFile -> idexo );
    msh.getElemNumMap          ( exoFile -> idexo );
    msh.getElementConnectivity ( exoFile -> idexo, mod );    
    // msh.getConnectivity        ( exoFile -> idexo );    
    ipl.findNodes              ( msh, mod, myfile );
    exo.writeNew               ( msh );
    
    // exit ( EXIT_SUCCESS );      
    
    // If the files containes refined nodes, estimate size.
    if ( msh.numFound != 0. )
    {                      
      estimateSize ( msh, size, willSplit, splitFile, exoFile -> fname );
      numRefine++;
    }
    
    // Destroy.
    msh.deallocateMesh     ( mod );
    exoFile -> closeFile   ( );       
  }    
  
  splitFile.close ();
  
  if ( willSplit == true )
  {
    std::cout << "\nHey bud, split your files before refining." 
      << " Do this by running ./scr/splitMesh.py" << std::endl;
  }
  else
  {
    std::cout << "\nLooks like the file doesn't need to be split!"
      << " Go ahead and run ./scr/refineElements why don't ya!" << std::endl;
  }
  
  return 0;
}

// FUNCTIONS.
void estimateSize ( Mesh &msh, int &size, bool &willSplit, 
  std::ofstream &splitFile, std::string fname )
{
  
  Constants con;
  
  double percent = (double (msh.numFound) / double (msh.num_nodes)) * 100.;
  std::cout << "Percentage of refined elements: " << percent << 
    "%" << std::endl;
  double newElem = msh.num_nodes - msh.numFound + msh.numFound * 
    pow ( 8, 2 );
  double newSize = ( (newElem+msh.num_elem) / msh.num_elem ) * size * con.toMB;

  std::cout << "Old size: "           << size * con.toMB << " MB" << std::endl;
  std::cout << "Estimated new size: " << newSize << " MB" << std::endl;
  
  if ( newSize > 250. )
  {
    std::cout << "I think " << fname << " will be larger than "
    << "250 MB after refinment. Splitting file ... " << std::endl;
    splitFile << fname << std::endl;
    willSplit = true;
  }
  
}

void setupRefineFile ( std::ofstream &myfile, std::string exoName, 
  int numRefine )
{
  
  std::string fRefineName = "./dat/refine/Refine";
  fRefineName.append ( std::to_string(numRefine) );
  fRefineName.append ( ".txt" );
  
  myfile.open ( fRefineName, std::ios::out );    
  myfile << exoName << std::endl;
  
}