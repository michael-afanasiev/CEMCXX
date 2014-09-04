#include <cmath>
#include "classes.hpp"
using namespace std;

void Utilities::pullInRad ( double &col, double &lon, double &rad, 
                 double &x,   double &y,   double &z, 
                 Mesh &msh )
{

  Utilities utl;

  if ( abs (rad - msh.radMax) < 1 )
    rad = rad - 1;

  if ( abs (rad - msh.radMin) < 1 )
    rad = rad + 1;

  utl.colLonRadRad2xyz ( col, lon, rad, x, y, z );

}

void Utilities::checkRegion ( Mesh &msh, double &rad )
{
  
  Constants con;
  
  if ( msh.regString == "innerCore" )
  {
    if ( rad >= con.innerCoreRad )
      rad = con.innerCoreRad - con.tiny;
  }

  if ( msh.regString == "outerCore" )
  {
    if ( rad >= (con.outerCoreRad-1) )
      rad = con.outerCoreRad - 1;
    if ( rad <= (con.innerCoreRad+1) )
      rad = con.innerCoreRad + 1;
  }
  
  if ( msh.regString == "lowerMantle" )
  {
    if ( rad >= con.R670 )
      rad = con.R670 - con.tiny;
    if ( rad <= con.outerCoreRad )
      rad = con.outerCoreRad + con.tiny;
  }
  
  if ( msh.regString == "transitionZone" )
  { 
    if (rad >= con.R400)
      rad = con.R400 - con.tiny;
    if (rad <= con.R670)
      rad = con.R670 + con.tiny;     
  }    
  
  if ( msh.regString == "upperMantle" )
  {
    if ( rad <= con.R400 )
      rad = con.R400 + con.tiny;
  }        
  
}

void Utilities::checkRegionExtr ( double x, double y, double z, short r,
                                  double &xUse, double &yUse, double &zUse )
{
  
  Constants con;
  Utilities utl;
  
  double col, lon, rad;
  
  utl.xyz2ColLonRadDeg ( x, y, z, col, lon, rad );
 
  /* Could change specfem to output three regions probably. Just like the normal
   * three regions. 'pullinrad' takes care of the rest i think. */

  if ( r == 1 )
  {
    if ( rad >= con.R_EARTH )
    {
      rad = con.R_EARTH - 1;
    }    
    if ( abs(rad-con.R020) <= con.tiny )
    {
      rad = con.R020 + 1;    
    }
    if ( abs(rad-con.R_EARTH) <= con.tiny )
    {
      rad = con.R_EARTH - 1;
    }
  }
   
  if ( r == 2 ) 
  {
    if ( abs(rad-con.R020) <= con.tiny )
    {
      rad = con.R020 - 1;
    }
    if ( abs(rad-con.R052) <= con.tiny )
    {
      rad = con.R052 + 1;
    }
  }
  
  if ( r == 3 ) 
  {
    if ( abs(rad-con.R052) <= con.tiny )
    {
      rad = con.R052 - 1;
    }
    if ( abs(rad-con.R100) <= con.tiny )
    {
      rad = con.R100 + 1;    
    }
  }
  
  if ( r == 4 ) 
  {
    if ( abs(rad-con.R100) <= con.tiny )
    {
      rad = con.R100 - 1;
    }
    if ( abs(rad-con.R400) <= con.tiny )
    {
      rad = con.R400 + 1;
    }
  }
  
  if ( r == 5 ) 
  {
    if ( abs(rad-con.R400) <= con.tiny )
    {
      rad = con.R400 - 1;
    }
    if ( abs(rad-con.R670) <= con.tiny )
    {
      rad = con.R670 + 1;
    }
  }
  
  if ( r == 6 ) 
  {
    if ( abs(rad-con.R670) <= con.tiny )
    {
      rad = con.R670 - 1;
    }
    if ( abs(rad-con.RTHO) <= con.tiny )
    {
      rad = con.RTHO + 1;
    }
  }

  if ( r == 7 ) 
  {
    if ( rad < con.outerCoreRad )
    {
      rad = con.outerCoreRad + 1;
    }
    
    if ( abs(rad-con.RTHO) <= con.tiny )
    {      
      rad = con.RTHO - 1;
    }
    if ( abs(rad-con.outerCoreRad) <= con.tiny )
    {
      rad = con.outerCoreRad + 1;
    }
  }
  
  if ( r == 8 ) 
  {
    if ( rad > con.outerCoreRad )
    {
      rad = con.outerCoreRad - 1;
    }
    if ( rad < con.innerCoreRad )
    {
      rad = con.innerCoreRad + 1;
    }
    if ( abs(rad-con.outerCoreRad) <= con.tiny )
    {
      rad = con.outerCoreRad - 1;
    }
    if ( abs(rad-con.innerCoreRad) <= con.tiny )
    {
      rad = con.innerCoreRad + 1;
    }
  }
  
  if ( r == 9 ) 
  {
    if ( rad > con.innerCoreRad )
    {
      rad = con.innerCoreRad - 5;
    }
    if ( rad < 0 )
    {
      rad = 5;
    }
    if ( abs(rad-con.innerCoreRad) <= con.bigtiny )
    {
      rad = con.innerCoreRad - 5;
    }
    if ( abs (rad) <= con.bigtiny )
    {
      rad = 5;      
    }
    
  }
  
  utl.colLonRadDeg2xyz ( col, lon, rad, xUse, yUse, zUse );
  
  
}

double Utilities::col2Lat ( double &col, char flag ) 
{
  /**
  Return the equivalent latitude from colatitude. The flag r or d describes 
  either radians or degrees.
  */
  
  Constants con;  
  double    lat=1e20;
  
  if ( flag == 'r' ) {
    lat = con.PIo2 - col;
  }
  else if ( flag == 'd' ) {
    lat = con.ninty - col;
  }
  
  if ( lat == 1e20 )
  {
    cout << "Something wrong in your lattitude conversion." << endl;
    exit ( EXIT_FAILURE );
  }
  
  return lat;
  
}

void Utilities::xyz2ColLonRadDeg ( double &x,   double &y,   double &z, 
                                   double &col, double &lon, double &rad )
{
  
  Constants con;
  
  if ( rad == 0 )
    rad = con.tiny;
  if ( z == 0 )
    z = con.tiny;
   
  rad = sqrt  ( x * x + y * y + z * z );
  col = acos  ( z / rad );
  lon = atan2 ( y, x );
  
  if ( rad > con.R_EARTH )
    rad = con.R_EARTH - 1;
  
  col = col * con.o80 / con.PI;
  lon = lon * con.o80 / con.PI;
  
}

void Utilities::xyz2ColLonRadRad ( double &x,   double &y,   double &z, 
                                   double &col, double &lon, double &rad )
{
  
  Constants con;
  
  if ( rad == 0 )
    rad = con.tiny;
  if ( z == 0 )
    z = con.tiny;
  
  rad = sqrt  ( x * x + y * y + z * z );
  col = acos  ( z / rad );
  lon = atan2 ( y, x );
  
  if ( rad > con.R_EARTH )
    rad = con.R_EARTH - 1;  
    
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
    mod.doRotate   = true;
    
  mod.colReg1    = false;
  mod.colReg2    = false;
  mod.lonReg1    = false;
  mod.lonReg2    = false;
  mod.lonReg3    = false;
  mod.lonReg4    = false;
  mod.wrapAround = false;      

  double a = mod.rotRad;
  double x = mod.rotVecX;
  double y = mod.rotVecY;
  double z = mod.rotVecZ;

  if ( mod.doRotate == true )
  {
    std::cout << "Rotation found. Rotating model " << mod.rotAng << 
      " degrees about ( " << x << ", " << y << ", " << z << " )." << std::endl;
  }
    
  double col, lon, rad;
  
  mod.forRot11 = cos(a) + (x * x) * (1 - cos(a));
  mod.forRot21 = z * sin(a) + x * y * (1 - cos(a));
  mod.forRot31 = y * sin(a) + x * z * (1 - cos(a));
  mod.forRot12 = x * y * (1 - cos(a)) - z * sin(a);
  mod.forRot22 = cos(a) + (y * y) * (1 - cos(a));
  mod.forRot32 = x * sin(a) + y * z * (1 - cos(a));
  mod.forRot13 = y * sin(a) + x * z * (1 - cos(a));
  mod.forRot23 = x * sin(a) + y * z * (1 - cos(a));
  mod.forRot33 = cos(a) + (z * x) * (1 - cos(a));
  
  mod.forRot23 = (-1) * mod.forRot23;
  mod.forRot31 = (-1) * mod.forRot31; 
  
  for ( size_t i=0; i<mod.x.size(); i++ ) 
  {

    double xNew;
    double yNew;
    double zNew;

    if( mod.doRotate == true ) {
    xNew = mod.forRot11 * mod.x[i] + mod.forRot21 * mod.y[i] + 
      mod.forRot31 * mod.z[i];
    yNew = mod.forRot12 * mod.x[i] + mod.forRot22 * mod.y[i] + 
      mod.forRot32 * mod.z[i];
    zNew = mod.forRot13 * mod.x[i] + mod.forRot23 * mod.y[i] + 
      mod.forRot33 * mod.z[i]; 
    } else {
      xNew = mod.x[i];
      yNew = mod.y[i];
      zNew = mod.z[i];
    }
    xyz2ColLonRadDeg ( xNew, yNew, zNew, col, lon, rad );
  
    if ( (col >=  0.) && (col <=  90.) )
      mod.colReg1 = true;
    if ( (col >= 90.) && (col <= 180.) )
      mod.colReg2 = true;
  
    if ( (lon >=    -1.) && (lon <=  91.) )
      mod.lonReg1 = true;      
    if ( (lon >=   89.) && (lon <= 181.) )
      mod.lonReg2 = true;      
    if ( (lon >=  179.) && (lon <= 271.) )
      mod.lonReg3 = true;      
    if ( (lon >=  269.) && (lon <= 361.) )
      mod.lonReg4 = true;      
    if ( (lon >=  359.) && (lon <= 451.) )
      mod.lonReg3 = true;
    if ( (lon >= -181.) && (lon <= -89.) )
      mod.lonReg3 = true;
    if ( (lon >=  -91.) && (lon <=   1.) )
      mod.lonReg4 = true;
    
    
    if ( rad <= 1221. && rad >= 0.    )
      mod.radReg1  = true;
    if ( rad <= 3480. && rad >= 1221. )
      mod.radReg2  = true;
    if ( rad <= 3571. && rad >= 3480. )
      mod.radReg3  = true;
    if ( rad <= 3671. && rad >= 3571. )
      mod.radReg4  = true;
    if ( rad <= 3771. && rad >= 3671. )
      mod.radReg5  = true;
    if ( rad <= 3871. && rad >= 3771. )
      mod.radReg6  = true;
    if ( rad <= 3971. && rad >= 3871. )
      mod.radReg7  = true;
    if ( rad <= 4071. && rad >= 3971. )
      mod.radReg8  = true;
    if ( rad <= 4171. && rad >= 4071. )
      mod.radReg9  = true;
    if ( rad <= 4271. && rad >= 4171. )
      mod.radReg10  = true;
    if ( rad <= 4371. && rad >= 4271. )
      mod.radReg11  = true;
    if ( rad <= 4471. && rad >= 4371. )
      mod.radReg12  = true;
    if ( rad <= 4571. && rad >= 4471. )
      mod.radReg13  = true;
    if ( rad <= 4671. && rad >= 4571. )
      mod.radReg14  = true;
    if ( rad <= 4771. && rad >= 4671. )
      mod.radReg15  = true;
    if ( rad <= 4871. && rad >= 4771. )
      mod.radReg16  = true;
    if ( rad <= 4971. && rad >= 4871. )
      mod.radReg17  = true;
    if ( rad <= 5071. && rad >= 4971. )
      mod.radReg18  = true;
    if ( rad <= 5171. && rad >= 5071. )
      mod.radReg19 = true;
    if ( rad <= 5271. && rad >= 5171. )
      mod.radReg20  = true;
    if ( rad <= 5371. && rad >= 5271. )
      mod.radReg21  = true;
    if ( rad <= 5426. && rad >= 5371. )
      mod.radReg22  = true;
    if ( rad <= 5481. && rad >= 5426. )
      mod.radReg23  = true;
    if ( rad <= 5536. && rad >= 5481. )
      mod.radReg24  = true;
    if ( rad <= 5591. && rad >= 5536. )
      mod.radReg25  = true;
    if ( rad <= 5646. && rad >= 5591. )
      mod.radReg26  = true;
    if ( rad <= 5701. && rad >= 5646. )
      mod.radReg27  = true;
    if ( rad <= 5746. && rad >= 5701. )
      mod.radReg28 = true;
    if ( rad <= 5791. && rad >= 5746. )
      mod.radReg29 = true;
    if ( rad <= 5836. && rad >= 5791. )
      mod.radReg30 = true;
    if ( rad <= 5881. && rad >= 5836. )
      mod.radReg31 = true;
    if ( rad <= 5926. && rad >= 5881. )
      mod.radReg32 = true;
    if ( rad <= 5971. && rad >= 5926. )
      mod.radReg33 = true;
    if ( rad <= 6021. && rad >= 5971. )
      mod.radReg34 = true;
    if ( rad <= 6071. && rad >= 6021. )
      mod.radReg35 = true;
    if ( rad <= 6121. && rad >= 6071. )
      mod.radReg36 = true;
    if ( rad <= 6171. && rad >= 6121. )
      mod.radReg37 = true;
    if ( rad <= 6221. && rad >= 6171. )
      mod.radReg38 = true;
    if ( rad <= 6271. && rad >= 6221. )
      mod.radReg39 = true;
    if ( rad <= 6291. && rad >= 6271. )
      mod.radReg40 = true;
    if ( rad <= 6311. && rad >= 6291. )
      mod.radReg41 = true;
    if ( rad <= 6331. && rad >= 6311. )
      mod.radReg42 = true;
    if ( rad <= 6351. && rad >= 6331. )
      mod.radReg43 = true;
    if ( rad <= 6361. && rad >= 6351. )
      mod.radReg44 = true;
    if ( rad <= 6371. && rad >= 6361. )
      mod.radReg45 = true;

  }
  
  /* Allow for the special case of wrapping along the -180/180 on the lon. axis.
  This is unnecessary if all 4 regions are selected. */
  
  if ( mod.lonReg2 == true && mod.lonReg3 == true )
  {
    mod.wrapAround = true;
  } else if  ( mod.lonReg1 == true && mod.lonReg2 == true 
            && mod.lonReg3 == true && mod.lonReg4 == true )
  {
    mod.wrapAround = false;
  }
    
  a             = ( -1 ) * mod.rotRad;
  mod.backRot11 = cos(a) + (x * x) * (1 - cos(a));
  mod.backRot21 = z * sin(a) + x * y * (1 - cos(a));
  mod.backRot31 = y * sin(a) + x * z * (1 - cos(a));
  mod.backRot12 = x * y * (1 - cos(a)) - z * sin(a);
  mod.backRot22 = cos(a) + (y * y) * (1 - cos(a));
  mod.backRot32 = x * sin(a) + y * z * (1 - cos(a));
  mod.backRot13 = y * sin(a) + x * z * (1 - cos(a));
  mod.backRot23 = x * sin(a) + y * z * (1 - cos(a));
  mod.backRot33 = cos(a) + (z * x) * (1 - cos(a));

  mod.backRot23 = (-1) * mod.backRot23;
  mod.backRot31 = (-1) * mod.backRot31;   
    
}                          

void Utilities::rotateForward ( double &xOld, double &yOld, double &zOld, 
                                double &xNew, double &yNew, double &zNew,
                                Model_file &mod )
{
  
  Constants con;
  
  if ( mod.doRotate == true ) 
  {  
    xNew = mod.forRot11 * xOld + mod.forRot21 * yOld + mod.forRot31 * zOld;
    yNew = mod.forRot12 * xOld + mod.forRot22 * yOld + mod.forRot32 * zOld;
    zNew = mod.forRot13 * xOld + mod.forRot23 * yOld + mod.forRot33 * zOld;    
  } 
  else 
  {  
    xNew = xOld;
    yNew = yOld;
    zNew = zOld; 
  }  
      
}

void Utilities::rotateBackward ( double &xOld, double &yOld, double &zOld, 
                                 double &xNew, double &yNew, double &zNew,
                                 Model_file &mod )
{
  
  Constants con;
  
  if ( mod.doRotate == true ) 
  {  
    xNew = mod.backRot11 * xOld + mod.backRot21 * yOld + mod.backRot31 * zOld;
    yNew = mod.backRot12 * xOld + mod.backRot22 * yOld + mod.backRot32 * zOld;
    zNew = mod.backRot13 * xOld + mod.backRot23 * yOld + mod.backRot33 * zOld;    
  } 
  else 
  {    
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

void Utilities::checkMeshEdge ( double &colOrig, double &lonOrig, Mesh &msh )

{
   
  if ( (colOrig < 1)    && (msh.colReg000_090 == true) )
    colOrig = 1;
  if ( (colOrig > 89 )  && (msh.colReg000_090 == true) )
    colOrig = 89;

//  if ( (colOrig < 91 )  && (msh.colReg090_180 == true) )
//    colOrig = 91;
  if ( (colOrig > 179)  && (msh.colReg090_180 == true) )
    colOrig = 179;

  if ( (lonOrig > 89)   && (msh.lonReg000_090 == true) )
    lonOrig = 89;
  if ( (lonOrig < 1 )   && (msh.lonReg000_090 == true) )
    lonOrig = 1;
  
  if ( (lonOrig < 91)   && (msh.lonReg090_180 == true) )
    lonOrig = 91;
  if ( (lonOrig > 179)  && (msh.lonReg090_180 == true) )
    lonOrig = 179;
  if ( (lonOrig < 0)    && (msh.lonReg090_180 == true) )
    lonOrig = 179;

  if ( (lonOrig > -1)   && (msh.lonReg270_360 == true) )
    lonOrig = -1;
  if ( (lonOrig < -89)  && (msh.lonReg270_360 == true) )
    
    lonOrig = -89;
  if ( (lonOrig > -91)  && (msh.lonReg180_270 == true) )
    lonOrig = -91;
  if ( (lonOrig < -179) && (msh.lonReg180_270 == true) )
    lonOrig = -179;
  if ( (lonOrig > 0)    && (msh.lonReg180_270 == true) )
    lonOrig = -179;
  
}
  
void Utilities::pullRad ( double &col, double &lon, double &rad, Mesh &msh, bool &fullSearch )
{

  bool     fixed = false;
  fullSearch     = false;
  double safeRad = 5.;
  double radPert = 0.1;

  if ( abs (rad - msh.radMax) < safeRad )
  {
    rad = rad - radPert;
  }

  if ( abs (rad - msh.radMax) == safeRad )
  {
    rad = msh.radMax - 2;
    fullSearch = true;
  }

  if ( (rad - msh.radMin) < safeRad )
  {
    rad = rad + radPert;
  }

  if ( abs (rad-msh.radMin) == safeRad )
  {
    rad=msh.radMin+2;
    fullSearch = true;
  }

}
  
  
