#include "classes.h"
#include <cmath>
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

void Model_file::populateSES3D ( string name, int &num_regions, 
  int &num_params, vector<vector<double>> &vec , char ftype ) {
    
  /** 
    
    This function will populate colatitude, longitude, and radius arrays 
    from an sed3d file. The values are stored in the model object.
    
  */
            
  string line;
  ifstream myfile;
  
  vector<double> sub;

  int * region;  
  int fix_stride;
  int fix_total;

  myfile.open ( name, ios::in );
  
  int lineno = 0;
  int inter  = 0;
  int regno  = 0;
  while ( getline (myfile, line) ) {
    
    // Read number of regions from first line in file.
    if ( lineno == 0 ) { 
      num_regions = stoi (line);
      region = new int [num_regions];
      
      // This is here to deal with the fact that the last coordinate in each
      // region is not assigned a parameter.
      if ( ftype == 'c' ) {
        fix_stride = 0;
        fix_total  = num_regions;
      } else if ( ftype == 'p' ) {
        fix_stride = 1;
        fix_total  = 0;
      }
      
    }
    
    // Read total number of parameters in file.
    if ( lineno == 1 ) {
      num_params = stoi (line);
    }
    
    // Read number of parameters in 1st region.
    if ( lineno == 2 ) {
      region[regno] = stoi (line);    
      regno++;
    }
    
    // Read rest of file.
    if ( lineno > 2 ) {
      inter++;
      
      // Stop just before the next region is encountered. If we are working
      // with a coordinate file, stop two before the next region is
      // encountered (there is no parameter associated with the last value).
      if ( inter < region[regno-1] + fix_stride ) {
        sub.push_back ( stod (line) );
      }
      
      // Suck up the number of parameters in the next region to continue.
      if ( inter == region[regno-1] + 1 ) {
        region[regno] = stoi (line);
                
        inter = 0;        
        regno++;
        
        vec.push_back ( sub );
        sub.clear ();        
      }
      
    }
     
    lineno++;
  }
  
  vec.push_back ( sub );
  sub.clear ();
  
  myfile.close ();
  
}

void Model_file::colLonRad2xyzSES3D ()
{
    
  int l = 0; 
  for ( int r=0; r<col_rad.size(); r++ ) {
    for ( int i=0; i<col_rad[r].size(); i++ ) {
      for ( int j=0; j<lon_rad[r].size(); j++ ) {
        for ( int k=0; k<rad[r].size();     k++ ) {
        
          x[l] = rad[r][k] * sin ( lon_rad[r][j] ) * cos ( col_rad[r][i] );
          y[l] = rad[r][k] * sin ( lon_rad[r][j] ) * sin ( col_rad[r][i] );
          z[l] = rad[r][k] * cos ( lon_rad[r][j] );
          l++;
        
        }
      }    
    }
  }
  
}

void Model_file::populateRadiansSES3D ()
{
  
  int i;
  
  vector<double> sub;
  
  Constants con;
  
  for ( int r=0; r<col_deg.size(); r++ ) {
    
    for ( i=0; i<col_deg[r].size(); i++ ) {
      sub.push_back ( col_deg[r][i] * con.PI / con.o80 );
    }
    
    col_rad.push_back ( sub );
    sub.clear ();
  
    for ( i=0; i<lon_deg[r].size(); i++ ) {
      sub.push_back ( lon_deg[r][i] * con.PI / con.o80 );;
    }
    
    lon_rad.push_back ( sub );
    sub.clear ();
    
  }
  
}

void Model_file::readSES3D ()
{
  
  int i;
  int n_head = 2;
  
  string imd = input_model_directory;
      
  // Generic -- we of course always need 3D coordinates.
  populateSES3D ( imd + "block_m_x", num_regions, num_x, col_deg, 'c' );
  populateSES3D ( imd + "block_m_y", num_regions, num_y, lon_deg, 'c' );
  populateSES3D ( imd + "block_m_z", num_regions, num_z, rad, 'c' );
  
  // Options for specific physics systems.
  if ( input_model_physics == "TTI" ) {
  populateSES3D ( imd + "dRHO", num_regions, num_p, rho, 'p' );
  populateSES3D ( imd + "dVSV", num_regions, num_p, vsv, 'p' );
  populateSES3D ( imd + "dVSH", num_regions, num_p, vsh, 'p' );
  populateSES3D ( imd + "dVPP", num_regions, num_p, vpp, 'p' );
  }
    
  // Put aside space for cartesian vectors (for speed).
  x.resize ( num_p );
  y.resize ( num_p );
  z.resize ( num_p );
    
}

void Model_file::readSPECFEM3D ()
{
  
  cout << "READING SPECFEM3D";
  
}