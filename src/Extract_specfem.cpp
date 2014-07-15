#include "classes.hpp"
#include <sstream>
#include <cstring>
#include <iomanip>
#include <string>

int main ( int argc, char *argv[] ) 
{

  bool dum;
  int regC  = atoi (argv[1]);
  int iProc = atoi (argv[2]);
  
  Region        reg;  
  Constants     con;
  Exodus_file   exo;
  Driver        drv;
  Model_file    mod;
  Utilities     utl;
  Discontinuity dis;
  
  std::cout << "Begin model building." << std::endl;
                        
  // Reset the random number generator (for use in the kdtree search)
  srand ( time (NULL) );      
  
  // Get the filename of the current cemRequest file.
  mod.getSpecFileName ( regC, iProc );
                
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
    msh.getInfo              ( exoFile -> idexo, 'p' );
    msh.getConnectivity      ( exoFile -> idexo );
    msh.createKDTreeUnpacked ( );

    if ( reg.colReg[fileIter] == "col000-090" )
      msh.colReg000_090 = true;

    if ( reg.colReg[fileIter] == "col090-180" )
      msh.colReg090_180 = true;

    if ( reg.lonReg[fileIter] == "lon000-090" )
      msh.lonReg000_090 = true;
     
    if ( reg.lonReg[fileIter] == "lon090-180" )
      msh.lonReg090_180 = true;

    if ( reg.lonReg[fileIter] == "lon180-270" )
      msh.lonReg180_270 = true;

    if ( reg.lonReg[fileIter] == "lon270-360" )
      msh.lonReg270_360 = true;

    fileIter += 1; 

    std::cout << "Extracting." << std::endl;
  
#pragma omp parallel for schedule (guided)        
    for ( size_t i=0; i<mod.x.size(); i++ )               
    {
      double c11, c12, c13, c14, c15, c16,c22, c23, c24, c25, c26, c33, c34;
      double c35, c36, c44, c45, c46, c55, c56, c66, rho, testX, testY, testZ;
      double col, lon, rad;  
      char skip;           
            
      utl.checkRegionExtr  ( mod.x[i], mod.y[i], mod.z[i], mod.r[i],
                             testX, testY, testZ ); 
    
                             
      utl.xyz2ColLonRadRad ( testX, testY, testZ, col, lon, rad ); 
      utl.fixTiny          ( testX, testY, testZ, col, lon, rad, skip, mod );  
      
      // TODO get rid of the stupid skip parameter.
      skip = 'p';


      if ( msh.lonMin < (-1 * con.PI / 2) && msh.lonMax > (con.PI / 2) )
	      msh.lonMax = -1 * con.PI / 2;

      if ( (rad <= msh.radMax) && 
           (rad >= msh.radMin) &&
           (lon <= (msh.lonMax + con.oneDegRad)) &&
           (lon >= (msh.lonMin - con.oneDegRad)) &&
           (col <= (msh.colMax + con.oneDegRad)) &&
           (col >= (msh.colMin - con.oneDegRad)) &&
            mod.r[i] != 9 ) 
      {            
    
        if ( abs (rad - msh.radMax) < 1 )
        {
            rad = rad - 1;
            utl.colLonRadRad2xyz ( col, lon, rad, testX, testY, testZ );
        }

        if ( abs (rad - msh.radMin) < 1 )
        {
            rad = rad + 1;
            utl.colLonRadRad2xyz ( col, lon, rad, testX, testY, testZ );
        }

        int pass = ipl.recover ( testX, testY, testZ, msh, c11, c12, c13, 
        c14, c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, 
        c45, c46, c55, c56, c66, rho, skip, dum ); 
          
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
      if ( mod.r[i] == 9 )
      {
        
        double vshUse, vppUse, rhoUse;
        
        Mod1d bm;
        bm.eumod ( rad, vshUse, vppUse, rhoUse);
          
        double vsvUse = vshUse;
        
        double N = rhoUse * vshUse * vshUse;
        double L = rhoUse * vsvUse * vsvUse;
        double A = rhoUse * vppUse * vppUse;
        
        double C = A;
        double F = A - 2 * L;
        double S = A - 2 * N;
        
        mod.c11[i] = C;
        mod.c22[i] = A;
        mod.c33[i] = A;
        mod.c12[i] = F;
        mod.c13[i] = F;
        mod.c23[i] = S;
        mod.c44[i] = N;
        mod.c55[i] = L;
        mod.c66[i] = L;
        mod.rhoUnwrap[i] = rhoUse;          
      }                          
    }        

    msh.deallocateMesh   ( mod );
    exoFile -> closeFile ( );                 
  }
  
  mod.projectSubspaceSPECFEM ( );

  std::cout << "Writing NetCDF" << std::endl;  
  mod.writeNetCDF (mod.rho, mod.specFileName+".rho");      
  mod.writeNetCDF (mod.vpp, mod.specFileName+".vpp");
  mod.writeNetCDF (mod.vsh, mod.specFileName+".vsh");
  mod.writeNetCDF (mod.vsv, mod.specFileName+".vsv");  

  return 0;
}
