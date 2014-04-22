#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <iomanip>

#include <netcdf>

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
  vector <vector <double> > &vec , char ftype ) {
    
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
    
    // Read number of parameters in 1st region.
    if ( lineno == 1 ) {
      region[regno] = stoi (line);    
      regno++;
    }
    
    // Read rest of file.
    if ( lineno > 1 ) {
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

void Model_file::createKDTreeUnpacked ( Mesh &msh )
{
  
  cout << "Creating KDTree ( model ).\n";
  tree  = kd_create (3);
  KDdat = new int [num_p];

#pragma omp parallel for
  for ( int i=0; i<num_p; i++ ) 
  { 
    KDdat[i] = i;
    if ( (radUnwrap[i] <= msh.radMax) && (radUnwrap[i] >= msh.radMin) )
    {   
      kd_insert3 ( tree, x[i], y[i], z[i], &KDdat[i] );
    }
  }
  
}

void Model_file::colLonRad2xyzSES3D ()
{
  
  Utilities util;
      
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
          
          radUnwrap[l] = rad[r][k];          
          
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
      
      if ( (lon_deg[r][i] >    0.) && (lon_deg[r][i] <  90.) )
        lonReg1 = true;      
      if ( (lon_deg[r][i] >   90.) && (lon_deg[r][i] < 180.) )
        lonReg2 = true;      
      if ( (lon_deg[r][i] >  180.) && (lon_deg[r][i] < 270.) )
        lonReg3 = true;      
      if ( (lon_deg[r][i] >  270.) && (lon_deg[r][i] < 360.) )
        lonReg4 = true;      
      if ( (lon_deg[r][i] >  360.) && (lon_deg[r][i] < 450.) )
        lonReg3 = true;
      if ( (lon_deg[r][i] > -180.) && (lon_deg[r][i] < -90.) )
        lonReg3 = true;
      if ( (lon_deg[r][i] >  -90.) && (lon_deg[r][i] <   0.) )
        lonReg4 = true;
      
      if ( lonReg2 == true && lonReg3 == true )
      {
        wrapAround = true;
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

void Model_file::populateRadians ( vector < vector <double> > &deg, 
                                   vector < vector <double> > &rad )
{
  
  Constants con;
  
  vector <double> sub;
  
  for ( int r=0; r<deg.size(); r++ ) {
    
    for ( int i=0; i<deg[r].size()-1; i++ ) {
      sub.push_back ( deg[r][i] * con.PI / con.o80 );
    }
    
    rad.push_back ( sub );
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
  
  if ( input_model_physics == "TTI" && intentions == "INTERPOLATE" ) {
    rhoUnwrap.resize ( num_p );
    vppUnwrap.resize ( num_p );
    vshUnwrap.resize ( num_p );
    vsvUnwrap.resize ( num_p );
  }
  
  if ( output_model_physics == "TTI" && intentions == "EXTRACT" ) {
    c11.resize ( num_p );
    c12.resize ( num_p );
    c13.resize ( num_p );
    c22.resize ( num_p );
    c23.resize ( num_p );
    c33.resize ( num_p );
    c44.resize ( num_p );
    c55.resize ( num_p );
    c66.resize ( num_p );    
    rhoUnwrap.resize ( num_p );
  }
    
  if ( intentions == "INTERPOLATE" ) {
    int l = 0;
    for ( int r=0; r<vsh.size(); r++ ) {
      for ( int i=0; i<vsh[r].size(); i++ ) {
      
        if ( input_model_physics == "TTI" ) {
                
          vsvUnwrap[l] = vsv[r][i];
          vshUnwrap[l] = vsh[r][i];
          vppUnwrap[l] = vpp[r][i];  
          rhoUnwrap[l] = rho[r][i];
                          
        }
      
        l++;
      }
    
    }
  }
  
}

void Model_file::projectSubspaceSPECFEM ( )
{
  
  vsh.resize(1);
  vsv.resize(1);
  vpp.resize(1);
  rho.resize(1);

  vsh[0].resize(c66.size());
  vsv[0].resize(c66.size());
  vpp[0].resize(c66.size());
  rho[0].resize(c66.size());
  
  int r = 0;
  for ( int i=0; i<c66.size(); i++ )
  {
    vsh[r][i] = sqrt ( c44[i] / rhoUnwrap[i] );
    vsv[r][i] = sqrt ( c55[i] / rhoUnwrap[i] );
    vpp[r][i] = sqrt ( c22[i] / rhoUnwrap[i] );
    rho[r][i] = rhoUnwrap[i];
  }
  
}

void Model_file::projectSubspaceSES3D ( )
{
  
  int ll = 0;
  for ( int r=0; r<col_deg.size(); r++ ) {
    
    int l  = 0;    
    for ( int i=0; i<col_rad[r].size(); i++ ) {
      for ( int j=0; j<lon_rad[r].size(); j++ ) {
        for ( int k=0; k<(rad[r].size()-1); k++ ) {
          
          if ( input_model_physics == "TTI" ) {
            vsh[r][l] = sqrt (c44[ll] / rhoUnwrap[ll]);
            vsv[r][l] = sqrt (c55[ll] / rhoUnwrap[ll]);
            vpp[r][l] = sqrt (c22[ll] / rhoUnwrap[ll]);
            rho[r][l] = rhoUnwrap[ll];
            l++;   
            ll++;     
          }
          
        }
      }
    }
  }            
  
}

void Model_file::readSES3D ()
{
  
  string imd = input_model_directory;
      
  // Generic -- we of course always need 3D coordinates.
  populateSES3D ( imd + "block_m_x", num_regions, col_deg, 'c' );
  populateSES3D ( imd + "block_m_y", num_regions, lon_deg, 'c' );
  populateSES3D ( imd + "block_m_z", num_regions, rad,     'c' );
  
  // Options for specific physics systems.  
  if ( intentions == "INTERPOLATE" ) {
    if ( input_model_physics == "TTI" ) {      
      populateSES3D ( imd + "dRHO", num_regions, rho, 'p' );
      populateSES3D ( imd + "dVSV", num_regions, vsv, 'p' );
      populateSES3D ( imd + "dVSH", num_regions, vsh, 'p' );
      populateSES3D ( imd + "dVPP", num_regions, vpp, 'p' );
    }
  } else if ( intentions == "EXTRACT" || intentions == "REFINE" ) {
    rho.resize ( num_regions );
    vsv.resize ( num_regions );
    vsh.resize ( num_regions );
    vpp.resize ( num_regions );       
  }
  
  // Calculate the proper number of parameters as there might be multiple
  // zones.
  num_p = 0;
  for ( int r=0; r!=num_regions; r++ ) 
  {  
    rho[r].resize ( (col_deg[r].size() - 1) * (lon_deg[r].size() - 1) * 
      (rad[r].size() - 1) );
    vsv[r].resize ( (col_deg[r].size() - 1) * (lon_deg[r].size() - 1) * 
      (rad[r].size() - 1) );
    vsh[r].resize ( (col_deg[r].size() - 1) * (lon_deg[r].size() - 1) * 
      (rad[r].size() - 1) );
    vpp[r].resize ( (col_deg[r].size() - 1) * (lon_deg[r].size() - 1) * 
      (rad[r].size() - 1) );         
    num_p += rho[r].size();
  }
  
  
  // Put aside space for cartesian vectors (for speed).
  x.resize ( num_p );
  y.resize ( num_p );
  z.resize ( num_p );
  radUnwrap.resize ( num_p );
    
  populateRadiansSES3D ();
  colLonRad2xyzSES3D   ();  
          
}

void Model_file::writeSES3D ()
{  
  string imd = input_model_directory;
  string omd = imd + "CEM";
  
  int status = mkdir ( omd.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
  
  if ( status != 0 )
  {
    cout << "Something fishy happened when I was creating the output " << 
      "directory. It's probably not a big deal ( i.e. the directory was " <<
      "just already present ) so I'm forging ahead." << endl;
  }
  
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
  
  using namespace netCDF;
  using namespace netCDF::exceptions;
  
  ifstream headerFile;
  string   line;
  
  int      lineNo;
  int      numCoord;
   
  try
  {
      
    NcFile dataFile ( specFileName, NcFile::read );
    
    NcVar dataX=dataFile.getVar("dataX");
    NcVar dataY=dataFile.getVar("dataY");
    NcVar dataZ=dataFile.getVar("dataZ"); 
    NcVar dataR=dataFile.getVar("regC_");       

    // Get the dim object so we can get the size of the data.
    NcDim dim = dataX.getDim (0);
    numCoord  = dim.getSize();
    
    double *dumX = new double [numCoord];
    double *dumY = new double [numCoord];
    double *dumZ = new double [numCoord];    
    short  *dumR = new short  [numCoord];
    
    num_p = numCoord;
    
    dataX.getVar (dumX);
    dataY.getVar (dumY);
    dataZ.getVar (dumZ);
    dataR.getVar (dumR);
    
    x.insert (x.end(), dumX, dumX+numCoord);
    y.insert (y.end(), dumY, dumY+numCoord);
    z.insert (z.end(), dumZ, dumZ+numCoord);
    r.insert (r.end(), dumR, dumR+numCoord);
    
    delete [] dumX;
    delete [] dumY;
    delete [] dumZ;
    delete [] dumR;
            
  } 
  catch (NcException& e)
  {
    
    e.what();
    cout << "FAILURE" << endl;
    
  }
  
  rho.resize ( 1 );
  vsv.resize ( 1 );
  vsh.resize ( 1 );
  vpp.resize ( 1 );   
  
  rho[0].resize ( numCoord );
  vsv[0].resize ( numCoord );
  vsh[0].resize ( numCoord );
  vpp[0].resize ( numCoord );
  
  rotAng = 0.;
  rotVecX = 0.;
  rotVecY = 0.;
  rotVecZ = 0.;      
}

int Model_file::writeNetCDF ( std::vector <std::vector<double>> &par, 
  std::string name )
{
  
  // Set up paramaeters for NetCDF writing ( taken from tutorial )
  static const int NC_ERR = 2;
  
  using namespace netCDF;
  using namespace netCDF::exceptions;
  
  // Compression filters ( 9 is highest )
  bool enableShuffleFilter = true;
  bool enableDeflateFilter = true;
  int deflateLevel         = 9;
  
  // Write to input / CEM directory.
  string imd = input_model_directory;
  string omd = imd + "CEM/";
  
  // Necessary try / catch block ( taken from tutorial )
  try
  {
    
    NcFile output ( name + ".nc", NcFile::replace );
        
    // Read total number of parameters.
    int totSize    = par[0].size();
                
    // Create new 1D array to output data from. Copy from packed regional arrays.
    double *dataOut = new double [totSize];    
    int l = 0;
    for ( int r=0; r<par.size(); r++ )
    {
      for ( vector <double> :: iterator p=par[r].begin(); p!=par[r].end(); 
        ++p )
      {
        // if ( *p < 1 )
        // {
        //   cout << "BAD " << l << endl;
        // }
        dataOut[l] = *p;
        l++;
      }
    }
    // Add a dimension to the netCDF file with n_par entries.
    NcDim dDim = output.addDim ( "param", totSize );
    
    // Create the parameter vector.
    vector <NcDim> dims;
    dims.push_back ( dDim );
    NcVar data = output.addVar ( "data", ncDouble, dims );
    
    // Set compression level.
    data.setCompression ( enableShuffleFilter, enableDeflateFilter, 
      deflateLevel );
      
    // Write to the data array.
    data.putVar ( dataOut );
    
    delete [] dataOut;
    
    return 0;    
    
  }
  catch ( NcException& e )
  {
    
    e.what ();
    return NC_ERR;
    
  }  
  
}

void Model_file::getSpecFileName ( int &regC, int &iProc )
{
  
  std::stringstream ssReg;
  std::stringstream ssPrc;
  
  ssReg << std::setw(2) << std::setfill('0');
  ssPrc << std::setw(4) << std::setfill('0');
        
  ssReg << std::to_string (regC+1);
  ssPrc << std::to_string (iProc);      
        
  specFileName = "./cemRequest/xyz_reg";
  specFileName.append (ssReg.str());
  specFileName.append ("_proc");
  specFileName.append (ssPrc.str());
  
}

void Model_file::deallocate ()
{
  
  kd_free ( tree );
  delete [] KDdat;
  
}
