#include "classes.hpp"
#include <sstream>
#include <cstring>
#include <iomanip>
#include <string>

int main ( int argc, char *argv[] ) 
{

  // Processor number of this particular file.
  int iProc = atoi (argv[1]);

  Region        reg;  
  Constants     con;
  Exodus_file   exo;
  Driver        drv;
  Model_file    mod;
  Utilities     utl;
  Discontinuity dis;
  
  std::cout << "Begin model building." << std::flush << std::endl;
                        
  // Reset the random number generator (for use in the kdtree search)
  srand ( time (NULL) );    

  // Get the filename of the current terragrid processor.
  mod.getTerraFileName ( iProc );

  std::cout << mod.terraFileName << std::endl;
  
  // Read parameter file.
  drv.initialize ( mod, dis, utl, exo, reg );    
  
  // Check usage
  drv.checkUsage ( mod, "EXTRACT" );

  std::cout << "\n----- Extracting -----\n";
 
  int fileIter = 0;
  for ( std::vector < Exodus_file > :: 
    iterator exoFile=reg.regionsExo.begin();
    exoFile!=reg.regionsExo.end(); ++exoFile ) 
  {
    std::cout << "\n";

    Mesh         msh;
    Interpolator ipl;

    exoFile -> openFile      ( exoFile -> fname );
    msh.getInfo              ( exoFile -> idexo );
    msh.getRegion            ( reg, fileIter );
    msh.getConnectivity      ( exoFile -> idexo );
    msh.createKDTreeUnpacked ( );

    std::cout << "Extracting." << std::flush << std::endl;
  
#pragma omp parallel for schedule (guided)        
    for ( size_t i=0; i<mod.x.size(); i++ )               
    {
      double c11, c12, c13, c14, c15, c16,c22, c23, c24, c25, c26, c33, c34;
      double c35, c36, c44, c45, c46, c55, c56, c66, rho, testX, testY, testZ;
      double col, lon, rad;  
            
      utl.checkRegionExtr  ( mod.x[i], mod.y[i], mod.z[i], mod.r[i],
                             testX, testY, testZ ); 
                                 
      utl.xyz2ColLonRadRad ( testX, testY, testZ, col, lon, rad ); 

      if ( (rad <= msh.radMax) && 
           (rad >= msh.radMin) &&
           (lon <= msh.lonMax) &&
           (lon >= msh.lonMin) &&
           (col <= msh.colMax) &&
           (col >= msh.colMin) )
             
      {            
    
        utl.pullInRad ( col, lon, rad, testX, testY, testZ, msh );
          
        int pass = ipl.recover ( testX, testY, testZ, msh, c11, c12, c13, 
        c14, c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, 
        c45, c46, c55, c56, c66, rho, 'p' ); 
          
        mod.c11[i]       = c11;
        mod.c12[i]       = c12;
        mod.c13[i]       = c13;
        mod.c22[i]       = c22;
        mod.c23[i]       = c23;
        mod.c33[i]       = c33;
        mod.c44[i]       = c44;
        mod.c55[i]       = c55;
        mod.c66[i]       = c66;            
        mod.rhoUnwrap[i] = rho;

      }
    }        

    msh.deallocateMesh   ( mod );
    exoFile -> closeFile ( );                 
    fileIter++; 
  }
  
  mod.projectSubspaceSPECFEM ( );
  mod.writeTERRAGRID ( );
  
  return 0;
}
