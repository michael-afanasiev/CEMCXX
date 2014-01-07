#include <iostream>
#include "classes.h"

void Exodus_file::openFile () {
    
  idexo = ex_open ( "/Users/michaelafanasiev/Development"
  "/SOURCE/CODE/comprehensive_earth_model/Exodus/single_scale_hex_highres.e", 
  EX_READ, &comp_ws, &io_ws, &vers ); 
  
  if (idexo < 0) {
    std::cout << "Fatal error opening exodus file. Exiting.\n";    
  } 
  
} 

void Exodus_file::closeFile () {
  
  ier = ex_close ( idexo );
  
  if (ier == 0) {
    std::cout << "File closed succesfully. \n";
  }
  
}