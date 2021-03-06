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
  
  if ( input_model_file_type == "SES3D" ) 
  {
    readSES3D ();    
  } 
  else if ( input_model_file_type == "SPECFEM3D" ) 
  {
    readSPECFEM3D ();    
  } else if ( input_model_file_type == "TERRAGRID" )
  {
    readTERRAGRID ();
  }
  else 
  {
    std::cout << "***MODEL FORMAT NOT RECOGNIZED. CHECK PARAM FILE. EXITING."
      << std::flush << std::endl;    
    exit (EXIT_FAILURE);  
  }
}

void Model_file::readTERRAGRID ()
{

  std::cout << "Reading TERRAGRID." << std::flush << std::endl;

  std::vector <double> terraX;
  std::vector <double> terraY;
  std::vector <double> terraZ;
  std::vector <double> terraRad;
  string line;
  ifstream radiusFile;
  ifstream xyzFile;

  string radName = "/mnt/lnec/afanasm/models/TERRAGRID/RADIUS";
  radiusFile.open     ( radName, ios::in );
  while ( getline (radiusFile, line) ) {

    terraRad.push_back ( stod (line) / 1000. );

  }

  radiusFile.close ();
  x.clear();
  y.clear();
  z.clear();

  double col1, col2, col3;
  string xyzName = terraFileName;
  xyzFile.open ( xyzName, ios::in );
  while ( getline (xyzFile, line) ) {

    std::istringstream ss(line);
    ss >> col1 >> col2 >> col3;
    terraX.push_back ( col1 );
    terraY.push_back ( col2 );
    terraZ.push_back ( col3 );

  }

  for ( size_t i=0; i<terraRad.size(); i++ ) {
    for ( size_t j=0; j<terraX.size(); j++ ) {

      x.push_back ( terraRad[i] * terraX[j] );
      y.push_back ( terraRad[i] * terraY[j] );
      z.push_back ( terraRad[i] * terraZ[j] );
      r.push_back ( 7 );

    }
  }


  rho.resize ( 1 );
  vsv.resize ( 1 );
  vsh.resize ( 1 );
  vpp.resize ( 1 );   
 
  num_p = terraRad.size()*terraX.size();

  rho[0].resize ( terraRad.size()*terraX.size() );
  vsv[0].resize ( terraRad.size()*terraX.size() );
  vsh[0].resize ( terraRad.size()*terraX.size() );
  vpp[0].resize ( terraRad.size()*terraX.size() );
  
  rotAng = 0.;
  rotVecX = 0.;
  rotVecY = 0.;
  rotVecZ = 0.;     

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

  int *region = 0;  
  int fix_stride;

  std::cout << "Reading: " << name << std::flush << std::endl;
    
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
    if ( lineno == 1 ) 
    {
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

  std::cout << "Creating KDTree ( model )." << std::flush << std::endl;

  std::vector <int> kdRegionsMin;
  std::vector <int> kdRegionsMax;

  /* Determine if we're crossing regions. */
  for ( size_t reg=0; reg<rad.size(); reg++ )
  {
    
    if ( ( msh.radMin <= maxRadReg[reg] && msh.radMin >= minRadReg[reg] ) )
      kdRegionsMin.push_back (reg);

    if ( ( msh.radMax <= maxRadReg[reg] && msh.radMax >= minRadReg[reg] ) )
      kdRegionsMax.push_back (reg);

  }

  bool diff = false;
  for ( size_t reg=0; reg<kdRegionsMin.size(); reg++ )
  {

    if ( kdRegionsMin[reg] != kdRegionsMax[reg] )
    {
      kdRegions.push_back ( kdRegionsMin[reg] );
      kdRegions.push_back ( kdRegionsMax[reg] );
      diff = true;
    } 

  }

  if      ( diff == false && kdRegionsMin.size() != 0 )
    kdRegions.push_back ( kdRegionsMin[0] );
  else if ( diff == false && kdRegionsMax.size() != 0 )
    kdRegions.push_back ( kdRegionsMax[0] );

  if ( kdRegions.size() == 0 )
  {
    std::cout << "__FATAL__ __FATAL__ __FATAL__  something bad in kdtree" 
      << std::flush <<std::endl;
    exit ( EXIT_FAILURE );
  }

  KDdat1 = new int [num_p];
  for ( int j=0; j<num_p; j++ )
    KDdat1[j] = j;
   
  KDdat2 = new int [num_p];
  for ( int j=0; j<num_p; j++ )
    KDdat2[j] = j;

  if ( kdRegions.size() <= 1 )
    tree1 = kd_create (3);
  if ( kdRegions.size() >= 2 )
    tree2 = kd_create (3);
  
  bool tree1Alive = false;
  bool tree2Alive = false;
  for ( int i=0; i<kdRegions.size(); i++ )
  {

#pragma omp parallel for
    for ( int j=0; j<num_p; j++ ) 
    { 

      if ( (radUnwrap[j] <= msh.radMax) && (radUnwrap[j] >= msh.radMin) && 
           (radUnwrap[j] <= maxRadReg[kdRegions[i]]) &&
           (radUnwrap[j] >= minRadReg[kdRegions[i]]) &&
           (i == 0) )
      {   
        kd_insert3 ( tree1, x[j], y[j], z[j], &KDdat1[j] );
        tree1Alive = true;
      }
      
      if ( (radUnwrap[j] <= msh.radMax) && (radUnwrap[j] >= msh.radMin) && 
           (radUnwrap[j] <= maxRadReg[kdRegions[i]]) &&
           (radUnwrap[j] >= minRadReg[kdRegions[i]]) &&
           (i == 1) )
      {   
        kd_insert3 ( tree2, x[j], y[j], z[j], &KDdat2[j] );
        tree2Alive = true;
      }
    }

  }

  if ( kdRegions.size () > 1 )
  {
    if ( tree1Alive == false )
    {
      kdRegions.erase (kdRegions.begin());
      tree1 = tree2;
    }
    if ( tree2Alive == false )
    {
      kdRegions.erase (kdRegions.begin()+1);
    }
  }

}


void Model_file::getMinMaxRegionSES3D ( )
{

  for ( size_t r=0; r <rad.size(); r++ )
  {
    minRadReg.push_back ( *min_element (rad[r].begin(), rad[r].end()) );
    maxRadReg.push_back ( *max_element (rad[r].begin(), rad[r].end()) );
  }

}

void Model_file::colLonRad2xyzSES3D ()
{
      
  /* IMPORTANT NOTE :: the range of k goes to rad.size-1, and col_rad and
  lon_rad don't, because col_rad and lon_rad have already been fixed up for
  their extra coordinate down below in populate radians. Rad hasn't, so it's
  done here. */
  int l = 0; 
  for ( size_t r=0; r<col_rad.size(); r++ ) {
    for ( size_t i=0; i<col_rad[r].size(); i++ ) {
      for ( size_t j=0; j<lon_rad[r].size(); j++ ) {
        for ( size_t k=0; k<rad[r].size()-1;     k++ ) {
       
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
  for ( size_t r=0; r<col_deg.size(); r++ ) 
  {
  
    for ( size_t i=0; i<col_deg[r].size()-1; i++ ) 
    {
      col_deg[r][i] = (col_deg[r][i] + col_deg[r][i+1]) / 2.;
      if ( col_deg[r][i] < colMin ) 
      {
        colMin = col_deg[r][i];
      } else if ( col_deg[r][i] > colMax ) 
      {
        colMax = col_deg[r][i];
      }
    }
    
    for ( size_t i=0; i<lon_deg[r].size()-1; i++ ) 
    {
      
      lon_deg[r][i] = (lon_deg[r][i] + lon_deg[r][i+1]) / 2.;


      // TODO __ MAKE THIS BETTER __
      if ( lon_deg[r][i] < lonMin ) 
      {
        lonMin = lon_deg[r][i];     
      } else if ( lon_deg[r][i] > lonMax ) 
      {
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
     
    for ( size_t i=0; i<rad[r].size()-1; i++ ) {
      rad[r][i] = (rad[r][i] + rad[r][i+1]) / 2.;
      if ( rad[r][i] < radMin ) {
        radMin = rad[r][i];
      } else if ( rad[r][i] > radMax ) {
        radMax = rad[r][i];
      }
      
    }
        
  }
    

  for ( size_t r=0; r<col_deg.size(); r++ ) {
    
    for ( size_t i=0; i<col_deg[r].size()-1; i++ ) {
      sub.push_back ( col_deg[r][i] * con.PI / con.o80 );
    }
    
    col_rad.push_back ( sub );
    sub.clear ();
    
  }
  
  for ( size_t r=0; r<lon_deg.size(); r++ ) {
  
    for ( size_t i=0; i<lon_deg[r].size()-1; i++ ) {
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
  
  for ( size_t r=0; r<deg.size(); r++ ) {
    
    for ( size_t i=0; i<deg[r].size()-1; i++ ) 
    {
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
    for ( size_t r=0; r<vsh.size(); r++ ) {
      for ( size_t i=0; i<vsh[r].size(); i++ ) {
      
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
  for ( size_t i=0; i<c66.size(); i++ )
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
  for ( size_t r=0; r<col_deg.size(); r++ ) {
    
    int l  = 0;    
    for ( size_t i=0; i<col_rad[r].size(); i++ ) {
      for ( size_t j=0; j<lon_rad[r].size(); j++ ) {
        for ( size_t k=0; k<(rad[r].size()-1); k++ ) {
          
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
  populateSES3D ( imd + "block_x", num_regions, col_deg, 'c' );
  populateSES3D ( imd + "block_y", num_regions, lon_deg, 'c' );
  populateSES3D ( imd + "block_z", num_regions, rad,     'c' );

  getMinMaxRegionSES3D ();
  

  // Options for specific physics systems.  
  if ( intentions == "INTERPOLATE" ) {
    
    if ( input_model_physics == "TTI" ) 
    {      
      populateSES3D ( imd + "dRHO", num_regions, rho, 'p' );
      populateSES3D ( imd + "dVSV", num_regions, vsv, 'p' );
      populateSES3D ( imd + "dVSH", num_regions, vsh, 'p' );
      populateSES3D ( imd + "dVPP", num_regions, vpp, 'p' );
    }

  } else if ( intentions == "EXTRACT" ) 
  {
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
  x.resize         ( num_p );
  y.resize         ( num_p );
  z.resize         ( num_p );
  radUnwrap.resize ( num_p );
    
  populateRadiansSES3D ( );
  colLonRad2xyzSES3D   ( );  
          
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

void Model_file::writeTERRAGRID ()
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
    dePopulateSES3D ( omd + "/rho." + terraFileProc, rho );
    dePopulateSES3D ( omd + "/vsv." + terraFileProc, vsv );
    dePopulateSES3D ( omd + "/vsh." + terraFileProc, vsh );
    dePopulateSES3D ( omd + "/vp."  + terraFileProc, vpp );    
  }
}

void Model_file::dePopulateSES3D ( string omf, vector < vector <double > > ou )
{
  
  std::cout << "Writing: " << omf << std::flush << std::endl;
  ofstream myfile ( omf, ios::out );
  
  myfile << ou.size() << "\n";
  for ( size_t r=0; r<ou.size(); r++ ) {
    myfile << ou[r].size() << "\n";
    for ( size_t i=0; i<ou[r].size(); i++ ) {
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
  
  int numCoord = 0;
   
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
    for ( size_t r=0; r<par.size(); r++ )
    {
      for ( vector <double> :: iterator p=par[r].begin(); p!=par[r].end(); 
        ++p )
      {
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

void Model_file::getTerraFileName ( int &iProc )
{

  std::stringstream ssPrc;

  ssPrc << std::setw(4) << std::setfill('0');
  ssPrc << std::to_string (static_cast<long long>(iProc));

  terraFileProc = ssPrc.str();
  terraFileName = "/mnt/lnec/afanasm/models/TERRAGRID/TerraGrid.";
  terraFileName.append (ssPrc.str());

}

void Model_file::getSpecFileName ( int &regC, int &iProc )
{
  
  std::stringstream ssReg;
  std::stringstream ssPrc;
  
  ssReg << std::setw(2) << std::setfill('0');
  ssPrc << std::setw(4) << std::setfill('0');
        
  ssReg << std::to_string (static_cast<long long>(regC+1));
  ssPrc << std::to_string (static_cast<long long>(iProc));      
        
  specFileName = "/mnt/lnec/afanasm/cemRequest/xyz_reg";
  specFileName.append (ssReg.str());
  specFileName.append ("_proc");
  specFileName.append (ssPrc.str());
  
}

void Model_file::deallocate ()
{

  delete [] KDdat1;
  delete [] KDdat2;
  kdRegions.clear();
  
}
