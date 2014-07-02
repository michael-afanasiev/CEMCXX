#include "classes.hpp"

void Mod1d::eumod ( double &rad, double &vsv, double &vpv, double &rho )
{
  
  Constants con;
  
  double x = rad / con.R_EARTH;

  /* Stretch mantle. */
  if ( (rad <= 6371) && (rad >= 6291) )  
  {
    rho = 2.6910 + 0.6924 * x;
    vpv = 4.1875 + 3.9382 * x - 0.035;
    vsv = 2.1519 + 2.3481 * x - 0.065;
  }
  else if ( (rad <= 6291) && (rad >= 6191) ) 
  {
    rho = 2.6910 + 0.6924 * x;
    vpv = 4.1875 + 3.9382 * x - 0.035;
    vsv = 2.1519 + 2.3481 * x - 0.065;
  }
  else if ( (rad <= 6191) && (rad >= 6051) ) 
  {
    rho = 9.1790 - 5.9841   * x;
    vpv = 40.5988 - 33.5317 * x - 0.035;
    vsv = 16.8261 - 12.7527 * x - 0.065;
  }
  else if ( (rad <= 6051) && (rad >= 5971) ) 
  {
    rho = 7.1089  - 3.8045  * x;
    vpv = 20.3926 - 12.2569 * x - 0.035;
    vsv = 8.9496  - 4.4597  * x - 0.065;
  }
  else if ( (rad <= 5971) && (rad >= 5771) )
  {
    rho = 11.2494 - 8.0298  * x;
    vpv = 39.7027 - 32.6166 * x - 0.035;
    vsv = 22.3512 - 18.5856 * x - 0.12;
  }
  else if ( (rad <= 5771) && (rad >= 5701) )   
  {
    rho = 5.3197  - 1.4836 * x;
    vpv = 19.0957 - 9.8672 * x - 0.035;
    vsv = 9.9839  - 4.9324 * x - 0.12;
  }
  else if ( (rad <= 5701) && (rad >= 5600) )   
  {
    rho = 7.9565  - 6.4761  * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 29.2766 - 23.6026 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 22.3459 - 17.2473 * x - 2.0834 * x * x + 0.9783 * x * x * x - 0.12;
  }
  else if ( (rad <= 5600) && (rad >= 3630) ) 
  {
    rho = 7.9565  - 6.4761  * x + 5.5283  * x * x -3.0807   * x * x * x;
    vpv = 24.9520 - 40.4673 * x + 51.4832 * x * x - 26.6419 * x * x * x;
    vsv = 11.1671 - 13.7818 * x + 17.4575 * x * x - 9.2777  * x * x * x - 0.12;
  }
  else if ( (rad <= 3630) && (rad >= 3480) )   
  {
    rho = 7.9565  - 6.4761 * x + 5.5283 * x * x - 3.0807 * x * x * x;
    vpv = 15.3891 - 5.3181 * x + 5.5242 * x * x - 2.5514 * x * x * x;
    vsv = 6.9254  + 1.4672 * x - 2.0834 * x * x + 0.9783 * x * x * x - 0.12;
  }
  else if ( (rad <= 3480) && (rad >= 1221.5) ) 
  {
     rho = 12.5815 - 1.2638 * x - 3.6426 * x * x - 5.5281  * x * x * x;
     vpv = 11.0487 - 4.0362 * x + 4.8023 * x * x - 13.5732 * x * x * x;
     vsv = 0.0;                                                        
  }
  else if (rad <= 1221.5) 
  {
     rho = 13.0885 - 8.8381 * x * x; 
     vpv = 11.2622 - 6.3640 * x * x; 
     vsv = 3.6678  - 4.4475 * x * x;
  }
  
}
