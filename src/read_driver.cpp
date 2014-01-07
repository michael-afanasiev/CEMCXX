#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
using namespace std;

void open_driver ( ifstream &myfile )
{
  
  myfile.open ("dat/driver.txt");
  return;
  
}

void get_token ( string &test )
{
  
  string del = "=";
  size_t pos = 0;
  
  pos  = test.find   (del);
  test = test.substr (pos + 2);
  return;
  
}

void read_driver_nMshFiles ( ifstream &myfile, string &num_mesh_files ) 
{  
  
  int i; 
  int ar_num;
  int LINENO_MSH_FILE = 2;  
  string line;

  i = 0;  
  if ( myfile.is_open() )
  {
    while ( getline (myfile, line) )
    {
      if ( i == LINENO_MSH_FILE ) 
      {
        get_token ( line );
      }
      i++;
    }
  }
  return;
    
}