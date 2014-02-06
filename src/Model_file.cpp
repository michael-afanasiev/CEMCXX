#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>

#include "classes.hpp"
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

  cout << "Reading: " << name << "\n";
    
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
        fix_stride = 1;
      } else if ( ftype == 'p' ) {
        fix_stride = 1;
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

void Model_file::populateParams ( Driver &drv, Exodus_file &exo ) 
{
  
  input_model_directory = drv.params[1];
  input_model_file_type = drv.params[2];
  input_model_physics   = drv.params[3];
  absolute_or_perturb   = drv.params[4]; 
  rotAng                = stod (drv.params[5]);
  rotVecX               = stod (drv.params[6]);
  rotVecY               = stod (drv.params[7]);
  rotVecZ               = stod (drv.params[8]);
  intentions            = drv.params[9];
  dimensions            = drv.params[10];
  interpolation         = drv.params[11];
  output_model_physics  = drv.params[12];
       
  cout << "\nModel information:\n* INPUT_MODEL_DIRECTORY: " <<
    input_model_directory << "\n* INPUT_MODEL_FILE_TYPE: "  <<
    input_model_file_type << "\n* INPUT_MODEL_PHYSICS: "    <<
    input_model_physics   << "\n* ABSOLUTE_OR_PERTURB: "    <<
    absolute_or_perturb   << "\n";
  
}

void Model_file::findMinMax ()
{
  
  double colMin = 180.0;
  double colMax = 0.0;
  
  for ( int r=0; r<col_deg.size(); r++ ) {
    for ( int i=0; i<col_deg[r].size(); i++ ) {
      if ( col_deg[r][i] < colMin ) {
        colMin = col_deg[r][i];
      }
      if ( col_deg[r][i] > colMax ) {
        colMax = col_deg[r][i];
      }
    }
  }
  
  double lonMin = 180.0;
  double lonMax = -180.0;
  
  for ( int r=0; r<lon_deg.size(); r++ ) {
    for ( int i=0; i<lon_deg[r].size(); i++ ) {
      if ( lon_deg[r][i] < lonMin ) {
        lonMin = lon_deg[r][i];
      }
      if ( lon_deg[r][i] > lonMax ) {
        lonMax = lon_deg[r][i];
      }
    }
  }
  
  double radMax = 0.0;
  double radMin = 6371.0;
  
  for ( int r=0; r<rad.size(); r++ ) {
    for ( int i=0; i<rad[r].size(); i++ ) {
      if ( rad[r][i] < radMin ) {
        radMin = rad[r][i];
      }
      if ( rad[r][i] > radMax ) {
        radMax = rad[r][i];
      }
    }
  }
  
  // cout << "Col: " << colMin << " " << colMax << "\n";
  // cout << "Lon: " << lonMin << " " << lonMax << "\n";
  // cout << "Rad: " << radMin << " " << radMax << "\n";
  
}

void Model_file::colLonRad2xyzSES3D ()
{
      
  /* IMPORTANT NOTE :: the range of k goes to rad.size-1, and col_rad and
  lon_rad don't, because col_rad and lon_rad have already been fixed up for
  their extra coordinate down below in populate radians. Rad hasn't, so it's
  done here. */
  int l = 0; 
  for ( int r=0; r<col_rad.size(); r++ ) {
    for ( int i=0; i<col_rad[r].size(); i++ ) {
      for ( int j=0; j<lon_rad[r].size(); j++ ) {
        for ( int k=0; k<rad[r].size()-1;     k++ ) {
        
          x[l] = rad[r][k] * cos ( lon_rad[r][j] ) * sin ( col_rad[r][i] );
          y[l] = rad[r][k] * sin ( lon_rad[r][j] ) * sin ( col_rad[r][i] );
          z[l] = rad[r][k] * cos ( col_rad[r][i] );
          l++;
        
        }
      }    
    }
  }
  
}

void Model_file::populateRadiansSES3D ()
{
  
  vector<double> sub;
  
  Constants con;
    
  // Compute box ceters and min/max on the fly.
  for ( int r=0; r<col_deg.size(); r++ ) {
  
    for ( int i=0; i<col_deg[r].size()-1; i++ ) {
      col_deg[r][i] = (col_deg[r][i] + col_deg[r][i+1]) / 2.;
      if ( col_deg[r][i] < colMin ) {
        colMin = col_deg[r][i];
      } else if ( col_deg[r][i] > colMax ) {
        colMax = col_deg[r][i];
      }
    }
    
    for ( int i=0; i<lon_deg[r].size()-1; i++ ) {
      lon_deg[r][i] = (lon_deg[r][i] + lon_deg[r][i+1]) / 2.;
      if ( lon_deg[r][i] < lonMin ) {
        lonMin = lon_deg[r][i];
      } else if ( lon_deg[r][i] > lonMax ) {
        lonMax = lon_deg[r][i];
      }
    }
    
    for ( int i=0; i<rad[r].size()-1; i++ ) {
      rad[r][i] = (rad[r][i] + rad[r][i+1]) / 2.;
      if ( rad[r][i] < radMin ) {
        radMin = rad[r][i];
      } else if ( rad[r][i] > radMax ) {
        radMax = rad[r][i];
      }
      
    }
        
  }

  for ( int r=0; r<col_deg.size(); r++ ) {
    
    for ( int i=0; i<col_deg[r].size()-1; i++ ) {
      sub.push_back ( col_deg[r][i] * con.PI / con.o80 );
    }
    
    col_rad.push_back ( sub );
    sub.clear ();
    
  }
  
  for ( int r=0; r<lon_deg.size(); r++ ) {
  
    for ( int i=0; i<lon_deg[r].size()-1; i++ ) {
      sub.push_back ( lon_deg[r][i] * con.PI / con.o80 );;
    }
    
    lon_rad.push_back ( sub );
    sub.clear ();
    
  }
  
}

void Model_file::openUp ( )
{
  
  /** This calculates the quantities of the elastic tensor that are not zero
  for a TTI model. 
  
  NOTE :: eta is set to 1 in this case. This may need to
  be changed in the future.
  */
  cout << "Changing physics.\n";
  
  if ( input_model_physics == "TTI" ) {
    c11    = new double [num_p]();
    c12    = new double [num_p]();
    c13    = new double [num_p]();
    c22    = new double [num_p]();
    c23    = new double [num_p]();
    c33    = new double [num_p]();
    c44    = new double [num_p]();
    c55    = new double [num_p]();
    c66    = new double [num_p]();  
    rhoMsh = new double [num_p]();
  }
    
  if ( intentions == "INTERPOLATE" ) {
    int l = 0;
    for ( int r=0; r<vsh.size(); r++ ) {
      for ( int i=0; i<vsh[r].size(); i++ ) {
      
        if ( input_model_physics == "TTI" ) {
                
          double N = rho[r][i] * vsh[r][i] * vsh[r][i];
          double L = rho[r][i] * vsv[r][i] * vsv[r][i];
          double A = rho[r][i] * vpp[r][i] * vpp[r][i];
    
          double C = A;
          double F = A - 2 * L;
      
          c11[l]    = C;
          c12[l]    = F;
          c13[l]    = F;
          c22[l]    = A;
          c23[l]    = A - 2 * N;
          c33[l]    = A;
          c44[l]    = N;
          c55[l]    = L;
          c66[l]    = L;            
          rhoMsh[l] = rho[r][i];
                
        }
      
        l++;
      }
    
    }
  }
  cout << "Physics has been changed.\n";  
  
}

void Model_file::readSES3D ()
{
  
  string imd = input_model_directory;
      
  // Generic -- we of course always need 3D coordinates.
  populateSES3D ( imd + "block_m_x", num_regions, num_x, col_deg, 'c' );
  populateSES3D ( imd + "block_m_y", num_regions, num_y, lon_deg, 'c' );
  populateSES3D ( imd + "block_m_z", num_regions, num_z, rad, 'c' );
  
  // Options for specific physics systems.  
  if ( intentions == "INTERPOLATE" ) {
    if ( input_model_physics == "TTI" ) {
    populateSES3D ( imd + "dRHO", num_regions, num_p, rho, 'p' );
    populateSES3D ( imd + "dVSV", num_regions, num_p, vsv, 'p' );
    populateSES3D ( imd + "dVSH", num_regions, num_p, vsh, 'p' );
    populateSES3D ( imd + "dVPP", num_regions, num_p, vpp, 'p' );
    }
  } else if ( intentions == "EXTRACT" ) {
    num_p = (num_x - 1) * (num_y - 1) * (num_z - 1);     
    rho.resize ( num_regions );
    vsv.resize ( num_regions );
    vsh.resize ( num_regions );
    vpp.resize ( num_regions );       
  }
  
  // Calculate the proper number of parameters as there might be multiple
  // zones.
  num_p = 0;
  for ( int r=0; r!=num_regions; r++ ) {
    
    rho[r].resize ( (col_deg[r].size() - 1) * (lon_deg[r].size() - 1) * 
      (rad[r].size() - 1) );
    vsv[r].resize ( (col_deg[r].size() - 1) * (lon_deg[r].size() - 1) * 
      (rad[r].size() - 1) );
    vsh[r].resize ( (col_deg[r].size() - 1) * (lon_deg[r].size() - 1) * 
      (rad[r].size() - 1) );
    vpp[r].resize ( (col_deg[r].size() - 1) * (lon_deg[r].size() - 1) * 
      (rad[r].size() - 1) );     
    num_p = rho[r].size() + num_p;
  }
  
  
  // Put aside space for cartesian vectors (for speed).
  x.resize ( num_p );
  y.resize ( num_p );
  z.resize ( num_p );
    
  populateRadiansSES3D ();
  colLonRad2xyzSES3D   ();  
      
      
}

void Model_file::writeSES3D ()
{  
  string imd = input_model_directory;
  string omd = imd + "CEM";
  
  int status = mkdir ( omd.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
  
  if ( output_model_physics == "TTI" ) {
    dePopulateSES3D ( omd + "/rho", rho );
    dePopulateSES3D ( omd + "/vsv", vsv );
    dePopulateSES3D ( omd + "/vsh", vsh );
    dePopulateSES3D ( omd + "/vp",  vpp );    
  }
}

void Model_file::dePopulateSES3D ( string omf, vector < vector <double > > ou )
{
  
  cout << "Writing: " << omf << "\n";
  ofstream myfile ( omf, ios::out );
  
  myfile << ou.size() << "\n";
  for ( int r=0; r<ou.size(); r++ ) {
    myfile << ou[r].size() << "\n";
    for ( int i=0; i<ou[r].size(); i++ ) {
      myfile << ou[r][i] << "\n";
    }
  }
  
  myfile.close();
  
}

void Model_file::readSPECFEM3D ()
{
  
  cout << "READING SPECFEM3D";
  
}