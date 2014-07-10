#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include "classes.hpp"
using namespace std;

void Driver::initialize ( Model_file &mod, Discontinuity &dis, Utilities &utl, 
  Exodus_file &exo, Region &reg )
{
  
  Constants con;
  
  int i;
  ifstream myfile;
  
  myfile.open ("./dat/driver.txt");
  
  if  (myfile.good () == !true ) {
    cout << "***Something is wrong with your param file. Exiting.\n";
    exit (EXIT_FAILURE);    
  }
  
  i = std::count ( istreambuf_iterator<char>(myfile), 
  istreambuf_iterator<char>(), '\n' );
    
  params = new string [i - N_HEADER];
  
  myfile.clear ();
  myfile.seekg (0, ios::beg);
  
  readDriver     ( myfile );
  closeDriver    ( myfile );
  populateParams ( mod );
    
  // Read model file.  
  if ( exo.allFiles == false )
    mod.read           ( );  
  
  // Read discontinuities.
  dis.read           ( );
  
  // Rotate model coordinates if necessary.
  utl.inquireRotate  ( mod );
      
  // Project to elastic tensor.
  if ( exo.allFiles == false )
    mod.openUp         ( );
  
  // Determine which exodus files to read.  
  exo.merge          ( reg, mod );
    
  return;
  
}

void Driver::populateParams ( Model_file &mod )
{
  
  Constants con;

  mod.mesh_directory        = params[0]; 
  mod.input_model_directory = params[1];
  mod.input_model_file_type = params[2];
  mod.input_model_physics   = params[3];
  mod.rotAng                = stod ( params[4] );
  mod.rotVecX               = stod ( params[5] );
  mod.rotVecY               = stod ( params[6] );
  mod.rotVecZ               = stod ( params[7] );
  mod.intentions            = params[8];
  mod.output_model_physics  = params[9];
  mod.subset                = params[10];
  
  if ( params[10] == "TRUE" )
  {
    mod.overwriteCrust = true;
  }
  else
  {
    mod.overwriteCrust = false;
  }

  if ( params[11] == "KERNEL1D" )
  {
    mod.kernel1d = true;
  }
  else if ( params[11] == "KERNEL3D" )
  {
    mod.kernel3d = true;
  }
  else
  {
    std::cout << "__OVERWRITING__" << std::flush << std::endl;
  }
  
  mod.rotRad = mod.rotAng * con.PI / con.o80;  
  
}

void Driver::closeDriver ( ifstream &myfile )
{
  
  myfile.close ();
  return;
  
}

void Driver::getToken ( string &test )
{
  
  string del = "=";
  size_t pos = 0;
  
  pos  = test.find   (del);
  test = test.substr (pos + 2);
  return;
  
}

void Driver::readDriver ( ifstream &myfile ) 
{  
  
  int i; 
  int j;
  string line;

  i = 0;  
  j = 0;
  if ( myfile.is_open() )
  {
    while ( getline (myfile, line) )
    {
      if ( i > N_HEADER ) 
      {
        getToken ( line );
        params[j] = line;
        j++;
      }
      i++;
    }
  }
  return;
    
}

void Driver::checkUsage ( Model_file &mod, std::string mode )
{
  
  if ( mod.intentions != mode ) 
  {
    
    std::cout << "\nHey bru I think you've got a mistake in your parameter " 
      << "file. I'm stopping here for you to check it out "
      << "( wrong intention? )" << std::endl;
    exit ( EXIT_FAILURE );
  }
  
}

void Driver::report ( Model_file &mod )
{
  
  ofstream myfile;
  
  time_t now = time (0);
  tm   *gmtm = gmtime (&now);
  char *dt   = asctime (gmtm);
  
  myfile.open ( "./dat/history/history.txt", std::ios::app );
  
  myfile << dt                        << std::endl;
  myfile << mod.input_model_directory << std::endl;
  myfile << mod.input_model_file_type << std::endl;
  myfile << mod.input_model_physics   << std::endl;
  myfile << mod.rotAng                << std::endl;
  myfile << mod.rotVecX               << std::endl;
  myfile << mod.rotVecY               << std::endl;
  myfile << mod.rotVecZ               << std::endl;    
  myfile << "\n"                      << std::endl;
  
  myfile.close ();

}
