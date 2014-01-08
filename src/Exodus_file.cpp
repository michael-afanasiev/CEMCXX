#include <stdio.h>
#include <string.h>
#include "classes.h"
using namespace std;

void Exodus_file::openFile ( string fname ) 
{
  
  char *cstr = new char [fname.length() + 1];
  strcpy ( cstr, fname.c_str() );
  
  cout << "Opening exodus file: " << cstr << "\n";
        
  idexo = ex_open ( cstr, EX_READ, &comp_ws, &io_ws, &vers ); 
  
  if (idexo < 0) {
    std::cout << "***Fatal error opening exodus file. Exiting.\n";  
    exit (EXIT_FAILURE);  
  } 
  
} 

void Exodus_file::closeFile () 
{
  
  ier = ex_close ( idexo );
  
  if (ier == 0) {
    std::cout << "File closed succesfully. \n";
  } else {
    cout << "***Fatal error closing exodus file. Exiting\n";
    exit (EXIT_FAILURE);  
  }
  
}