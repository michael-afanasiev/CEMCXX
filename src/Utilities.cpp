#include <cmath>
#include "classes.hpp"
using namespace std;

double Utilities::col2Lat ( double &col, char flag ) 
{
  /**
  Return the equivalent latitude from colatitude. The flag r or d describes 
  either radians or degrees.
  */
  
  Constants con;  
  double    lat;
  
  if ( flag == 'r' ) {
    lat = con.PIo2 - col;
  }
  else if ( flag == 'd' ) {
    lat = con.ninty - col;
  }
  
  return lat;
  
}

void Utilities::xyz2ColLonRadDeg ( double &x,   double &y,   double &z, 
                                   double &col, double &lon, double &rad )
{
  
  Constants con;
  
  rad = sqrt  ( x * x + y * y + z * z );
  col = acos  ( z / rad ) * con.o80 / con.PI ;
  lon = atan2 ( y, x ) * con.o80 / con.PI ;
  
}

void Utilities::xyz2ColLonRadRad ( double &x,   double &y,   double &z, 
                                   double &col, double &lon, double &rad )
{
  
  Constants con;
  
  rad = sqrt  ( x * x + y * y + z * z );
  col = acos  ( z / rad );
  lon = atan2 ( y, x );
  
}
                                  
void Utilities::colLonRadDeg2xyz ( double col,  double lon,  double rad,
                                   double &x,   double &y,   double &z ) 
{
  
  Constants con;
    
  
  col = col * con.PI / con.o80;
  lon = lon * con.PI / con.o80;
  
  x = rad * cos (lon) * sin (col);
  y = rad * sin (lon) * sin (col);
  z = rad * cos (col);  
  
}    

void Utilities::colLonRadRad2xyz ( double col,  double lon,  double rad,
                                   double &x,   double &y,   double &z ) 
{
  
  Constants con;
  
  x = rad * cos (lon) * sin (col);
  y = rad * sin (lon) * sin (col);
  z = rad * cos (col);  
  
}                                

void Utilities::rotate ( Model_file &mod )
{
  /**
  Just a classic rotation matrix. Give me x, y, z and I'll give x_out, y_out
  z_out
  */
  
  Constants con;
  
  if ( mod.rotAng != 0 )
  cout << "Rotating model by " << mod.rotAng << " degrees about (" <<
    mod.rotVecX << ", " << mod.rotVecY << ", " << mod.rotVecZ << ")\n";
  
  double a = mod.rotRad;
  double x = mod.rotVecX;
  double y = mod.rotVecY;
  double z = mod.rotVecZ;
    
  double rot11 = cos(a) + (x * x) * (1 - cos(a));
  double rot22 = cos(a) + (y * y) * (1 - cos(a));
  double rot33 = cos(a) + (z * z) * (1 - cos(a));  
  double rot21 = z * sin(a) + x * y * (1 - cos(a));
  double rot31 = y * sin(a) + x * z * (1 - cos(a));
  double rot32 = x * sin(a) + y * z * (1 - cos(a));
  double rot13 = y * sin(a) + x * z * (1 - cos(a));
  double rot23 = x * sin(a) + y * z * (1 - cos(a));
  double rot12 = x * y * (1 - cos(a)) - z * sin(a);
  
  mod.colMin = 180.;
  mod.colMax = 0.;
  mod.lonMin = 180.;
  mod.lonMax = -180.;
  mod.radMin = 6371.;
  mod.radMax = 0.;
  
  
  for ( int i=0; i<mod.x.size (); i++ ) {  
    mod.x[i] = rot11 * mod.x[i] + rot12 * mod.y[i] + rot13 * mod.z[i];
    mod.y[i] = rot21 * mod.x[i] + rot22 * mod.y[i] + rot23 * mod.z[i];
    mod.z[i] = rot31 * mod.x[i] + rot32 * mod.y[i] + rot33 * mod.z[i]; 
    
    double col, lon, rad;  
    
    xyz2ColLonRadDeg ( mod.x[i], mod.y[i], mod.z[i], col, lon, rad );
    
    if ( (mod.x[i] != 0) || (mod.y[i] != 0) ) {

      if ( col < mod.colMin )
        mod.colMin = col;
      if ( col > mod.colMax )
        mod.colMax = col;
      if ( lon < mod.lonMin )
        mod.lonMin = lon;
      if ( lon > mod.lonMax )
        mod.lonMax = lon;
      if ( rad < mod.radMin )
        mod.radMin = rad;
      if ( rad > mod.radMax )
        mod.radMax = rad;
    
    }

    if ( (mod.x[i] == 0) && (mod.y[i] == 0) && (mod.z[i] > 0) )
      mod.colMin = 0;
    if ( (mod.x[i] == 0) && (mod.y[i] == 0) && (mod.z[i] < 0) )
      mod.colMax = 180;
    
  }
      
}

void Utilities::convertBary ( double &xp, double &yp, double &zp, 
                              double &x1, double &x2, double &x3, double &x4,
                              double &y1, double &y2, double &y3, double &y4,
                              double &z1, double &z2, double &z3, double &z4,
                              double &l1, double &l2, double &l3, double &l4 ) 
{
  
  double vecX = xp - x4;
  double vecY = yp - y4;
  double vecZ = zp - z4;
  
  double a = x1 - x4;
  double d = y1 - y4;
  double g = z1 - z4;
  double b = x2 - x4;
  double e = y2 - y4;
  double h = z2 - z4;
  double c = x3 - x4;
  double f = y3 - y4;
  double i = z3 - z4;
  
  double det = 1 / (( a * ( e * i - f * h ) ) - ( b * ( i * d - f * g ) ) +
    ( c * ( d * h - e * g ) ));
  
  double ai = det * (e * i - f * h);
  double bi = det * (d * i - f * g) * (-1);
  double ci = det * (d * h - e * g);
  double di = det * (b * i - c * h) * (-1);
  double ei = det * (a * i - c * g);
  double fi = det * (a * h - b * g) * (-1);
  double gi = det * (b * f - c * e);
  double hi = det * (a * f - c * d) * (-1);
  double ii = det * (a * e - b * d);
  
  l1 = ai * vecX + di * vecY + gi * vecZ;
  l2 = bi * vecX + ei * vecY + hi * vecZ;
  l3 = ci * vecX + fi * vecY + ii * vecZ;
  l4 = 1 - l1 - l2 - l3;
  
}