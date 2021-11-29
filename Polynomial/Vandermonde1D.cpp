// Vandermonde1D.m
// function [V1D] = Vandermonde1D(N,xp)
// 2007/06/06
//---------------------------------------------------------
#include "Jacobi1D.h"

//---------------------------------------------------------
void Vandermonde1D(const int &N, const VectorXd& r,
                   MatrixXd& V1D)
//---------------------------------------------------------
{
  // function [V1D] = Vandermonde1D(N,xp)
  // Purpose : Initialize the 1D Vandermonde Matrix.
  //	    V_{ij} = phi_j(xp_i);

  int Nc = r.size();
  VectorXd P;
  V1D = MatrixXd::Constant(Nc, N+1, 0.);
  for (int j=0; j<=N; ++j) {
     JacobiP(r, 0, 0, j, P);
      V1D.col(j) = P;
  }
  return ;
}
