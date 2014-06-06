#include "classes.hpp"
  
void Tensor::TTI ( Mesh &msh, Model_file &mod, 
                   double &N, double &L, double &A, double &C, double &F, 
                   double &S, double &rhoUse, double &qvCor, double &tap )
{
  
  Mod1d bm;
  
  // Get 1d background values.
  double vs1d, vp1d, rho1d;
  bm.eumod                   ( mshRadRot, vs1d, vp1d, rho1d );
  
  double vshExo = sqrt ( msh.c44[i] / msh.rho[i] );
  double vsvExo = sqrt ( msh.c55[i] / msh.rho[i] );
  double vppExo = sqrt ( msh.c22[i] / msh.rho[i] );
  double rhoExo = msh.rho[i];
              
  double vshMod = mod.vshUnwrap[point];
  double vsvMod = mod.vsvUnwrap[point];
  double vppMod = mod.vppUnwrap[point];
  double rhoMod = mod.rhoUnwrap[point]; 

  double rhoModCor;
  double vshModCor;
  double vsvModCor;
  double vppModCor;

  if ( mod.kernel == false )
  {        
    rhoModCor = rhoMod;
    vshModCor = vshMod * qvCor;
    vsvModCor = vsvMod * qvCor;
    vppModCor = vppMod * qvCor;
  }
  else
  {
    rhoModCor = ( rhoMod + rho1d );      
    vshModCor = ( vshMod + vs1d ) * qvCor;
    vsvModCor = ( vsvMod + vs1d ) * qvCor;
    vppModCor = ( vppMod + vp1d ) * qvCor;
  }
        
  double rhoUse = ( 1 - tap ) * rhoExo + tap * rhoModCor;
  double vshUse = ( 1 - tap ) * vshExo + tap * vshModCor;
  double vsvUse = ( 1 - tap ) * vsvExo + tap * vsvModCor;
  double vppUse = ( 1 - tap ) * vppExo + tap * vppModCor;
        
  N = rhoUse * vshUse * vshUse;
  L = rhoUse * vsvUse * vsvUse;
  A = rhoUse * vppUse * vppUse;
  
  C = A;
  F = A - 2 * L;
  S = A - 2 * N;
  
}