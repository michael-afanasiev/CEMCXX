#include "classes.h"
using namespace std;

void Model_file::read () 
{
  
  if ( input_model_file_type == "SES3D" ) {
    readSES3D ();    
  } else if ( input_model_file_type == "SPECFEM3D") {
    readSPECFEM3D ();    
  } else {
    cout << "***MODEL FORMAT NOT RECOGNIZED. CHECK PARAM FILE. EXITING.\n";
    exit (EXIT_FAILURE);  
  }
  
}

void populateSES3D ( string name, int &num_regions, int &num_params, 
  vector<double> &vec , char ftype ) {
            
  string line;
  ifstream myfile;

  int * region;  
  int fix;

  // This is here to deal with the fact that the last coordinate in each
  // region is not assigned a parameter.
  if ( ftype == 'c' ) {
    fix = 0;
  } else if ( ftype == 'p' ) {
    fix = 1;
  }

  myfile.open ( name );
  
  int lineno = 0;
  int inter  = 0;
  int regno  = 0;
  while ( getline (myfile, line) ) {
    
    if ( lineno == 0 ) { 
      num_regions = stoi (line);
      region = new int [num_regions];
    }
    
    if ( lineno == 1 ) {
      num_params = stoi (line);
    }
    
    if ( lineno == 2 ) {
      region[regno] = stoi (line);    
      regno++;
    }
    
    if ( lineno > 2 ) {
      inter++;
      
      if ( inter < region[regno-1] + fix ) {
        vec.push_back (stod (line) );
      }
      
      if ( inter == region[regno-1] + 1 ) {
        region[regno] = stoi (line);
        inter = 0;
        regno++;
      }
      
    }
     
    lineno++;
  }
  
}

void Model_file::readSES3D ()
{
  
  int i;
  int n_head = 2;
  
  string imd = input_model_directory;
      
  cout << "READING SES3D.";
  populateSES3D ( imd + "block_m_x", num_regions, num_x, x, 'c' );
  populateSES3D ( imd + "block_m_y", num_regions, num_y, y, 'c' );
  populateSES3D ( imd + "block_m_z", num_regions, num_z, z, 'c' );
  
  if ( input_model_physics == "TTI" ) {
  populateSES3D ( imd + "dRHO", num_regions, num_p, rho, 'p' );
  populateSES3D ( imd + "dVSV", num_regions, num_p, vsv, 'p' );
  populateSES3D ( imd + "dVSH", num_regions, num_p, vsh, 'p' );
  populateSES3D ( imd + "dVPP", num_regions, num_p, vpp, 'p' );
  }
  
}

void Model_file::readSPECFEM3D ()
{
  
  cout << "READING SPECFEM3D";
  
}