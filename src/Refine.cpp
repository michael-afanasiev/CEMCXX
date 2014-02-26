#include "classes.hpp"
#include "kdtree.h"

int main () {

  std::cout << "Begin refinement procedure.\n";

  std::ifstream myfile;

  Constants    con;
  Exodus_file  exo;
  Driver       drv;
  Mesh         msh;
  Model_file   mod;
  Utilities    util;
  Interpolator inter;
  
  double c11, c12, c13, c14, c15, c16,c22, c23, c24, c25, c26, c33, c34, c35;
  double c36, c44, c45, c46, c55, c56, c66, rho, testX, testY, testZ;

  std::cout << "Begin model building.\n";  

  // ********************************************************************* //
  //                       READING THE PARAMETER FILE                      //
  // ********************************************************************* //
  
  // Determine how many lines to read from parameter file.
  drv.openDriver ( myfile );

  // Read those lines from parameter file and put them into 'params'.
  drv.readDriver ( myfile );

  // Close driver.
  drv.closeDriver ( myfile );

  // ********************************************************************* //
  //                            READ IN PARAMS                             //
  // ********************************************************************* //

  mod.populateParams ( drv, exo );      

  // ********************************************************************* //
  //                            READ IN MODEL                              //
  // ********************************************************************* //


  mod.read ();  
  
  exo.openFile ( "./dat/test.ex2" );
  
  msh.getInfo ( exo.idexo, 'c' );    
  msh.getConnectivity ( exo.idexo );
  
  std::cout << "Creating KDTree.\n";
  kdtree *tree = kd_create (3);
  int *dat     = new int [msh.num_nodes];
  for ( int i=0; i<msh.num_nodes; i++ ) {
    dat[i] = i;
    kd_insert3 ( tree, msh.xmsh[i], msh.ymsh[i], msh.zmsh[i], &dat[i] );
  }    
  
  std::cout << "Done.\n";
  
  
  std::ofstream outfile ("elements.txt");
  for ( int r=0; r<mod.col_rad.size(); r++ ) {
    for ( int i=0; i<mod.col_rad[r].size(); i++ ) {
      std::cout << "Col number: " << i << "\n";    
      
      for ( int j=0; j<mod.lon_rad[r].size(); j++ ) {
        for ( int k=0; k<mod.rad[r].size()-1;     k++ ) {
                    
          util.colLonRadRad2xyz ( mod.col_rad[r][i], mod.lon_rad[r][j], 
            mod.rad[r][k], testX, testY, testZ );
                        
          int element = inter.recover ( testX, testY, testZ, tree, msh,
              c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26,
              c33, c34, c35, c36, c44, c45, c46, c55, c56, c66,
              rho, 'e' );
              
          outfile << element << "\n";
        }
      }
      
    }
  }
   
}