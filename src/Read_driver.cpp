#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include "classes.hpp"
using namespace std;

void Driver::openDriver ( ifstream &myfile )
{
  
  int i;
  
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
  
  readDriver  ( myfile );
  closeDriver ( myfile );
  
  return;
  
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