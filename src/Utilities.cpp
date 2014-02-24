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

void Utilities::inquireRotate ( Model_file &mod )
{
  
  if ( mod.rotAng != 0 ) 
  {
    
    mod.doRotate = true;
    
    double a = mod.rotRad;
    double x = mod.rotVecX;
    double y = mod.rotVecY;
    double z = mod.rotVecZ;

    std::cout << "Rotation found. Rotating model " << mod.rotAng << 
      " degrees about ( " << x << ", " << y << ", " << z << " )." << std::endl;
  
    double colMin = 180.;
    double colMax = 0.;
    double lonMin = 360.;
    double lonMax = -360.;
    double radMin = 6371.;
    double radMax = 0.;
  
    double col, lon, rad;
  
    mod.rot11 = cos(a) + (x * x) * (1 - cos(a));
    mod.rot21 = z * sin(a) + x * y * (1 - cos(a));
    mod.rot31 = y * sin(a) + x * z * (1 - cos(a));
    mod.rot12 = x * y * (1 - cos(a)) - z * sin(a);
    mod.rot22 = cos(a) + (y * y) * (1 - cos(a));
    mod.rot32 = x * sin(a) + y * z * (1 - cos(a));
    mod.rot13 = y * sin(a) + x * z * (1 - cos(a));
    mod.rot23 = x * sin(a) + y * z * (1 - cos(a));
    mod.rot33 = cos(a) + (z * x) * (1 - cos(a));
  
    mod.rot23 = (-1) * mod.rot23;
    mod.rot31 = (-1) * mod.rot31; 
  
    for ( int i=0; i<mod.x.size(); i++ ) 
    {
    
      double xNew = mod.rot11 * mod.x[i] + mod.rot21 * mod.y[i] + 
        mod.rot31 * mod.z[i];
      double yNew = mod.rot12 * mod.x[i] + mod.rot22 * mod.y[i] + 
        mod.rot32 * mod.z[i];
      double zNew = mod.rot13 * mod.x[i] + mod.rot23 * mod.y[i] + 
        mod.rot33 * mod.z[i]; 
    
      xyz2ColLonRadDeg ( xNew, yNew, zNew, col, lon, rad );
      
      if ( (col >  0.) && (col <  90.) )
        mod.colReg1 = true;
      if ( (col > 90.) && (col < 180.) )
        mod.colReg2 = true;
      
      if ( (lon >    0.) && (lon <  90.) )
        mod.lonReg1 = true;      
      if ( (lon >   90.) && (lon < 180.) )
        mod.lonReg2 = true;      
      if ( (lon >  180.) && (lon < 270.) )
        mod.lonReg3 = true;      
      if ( (lon >  270.) && (lon < 360.) )
        mod.lonReg4 = true;      
      if ( (lon >  360.) && (lon < 450.) )
        mod.lonReg3 = true;
      if ( (lon > -180.) && (lon < -90.) )
        mod.lonReg3 = true;
      if ( (lon >  -90.) && (lon <   0.) )
        mod.lonReg4 = true;
      
    }
        
    a         = ( -1 ) * mod.rotRad;
    mod.rot11 = cos(a) + (x * x) * (1 - cos(a));
    mod.rot21 = z * sin(a) + x * y * (1 - cos(a));
    mod.rot31 = y * sin(a) + x * z * (1 - cos(a));
    mod.rot12 = x * y * (1 - cos(a)) - z * sin(a);
    mod.rot22 = cos(a) + (y * y) * (1 - cos(a));
    mod.rot32 = x * sin(a) + y * z * (1 - cos(a));
    mod.rot13 = y * sin(a) + x * z * (1 - cos(a));
    mod.rot23 = x * sin(a) + y * z * (1 - cos(a));
    mod.rot33 = cos(a) + (z * x) * (1 - cos(a));

    mod.rot23 = (-1) * mod.rot23;
    mod.rot31 = (-1) * mod.rot31;   
    
  }
    
}                          

void Utilities::rotate ( double &xOld, double &yOld, double &zOld, 
                         double &xNew, double &yNew, double &zNew,
                         Model_file &mod )
{
  
  Constants con;
  
  if ( mod.doRotate == true ) {
    
    xNew = mod.rot11 * xOld + mod.rot21 * yOld + mod.rot31 * zOld;
    yNew = mod.rot12 * xOld + mod.rot22 * yOld + mod.rot32 * zOld;
    zNew = mod.rot13 * xOld + mod.rot23 * yOld + mod.rot33 * zOld; 
      
  } else {
    
    xNew = xOld;
    yNew = yOld;
    zNew = zOld;
    
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