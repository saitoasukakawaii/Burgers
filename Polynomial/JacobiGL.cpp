//---------------------------------------------------------
#include "Jacobi1D.h"
#include "cassert"
//---------------------------------------------------------
void JacobiGL(const double &alpha, const double &beta, const int &N,
                VectorXd& x)
//---------------------------------------------------------
{
  // function [x] = JacobiGL(alpha,beta,N)
  // Purpose: Compute the N'th order Gauss Lobatto quadrature 
  //          points, x, associated with the Jacobi polynomial,
  //          of type (alpha,beta) > -1 ( <> -0.5). 
  // N: order of polynomial, return N+1 points of Gauss Lobatto quadrature
  assert(N>0);
  x = VectorXd::Constant(N+1, 0.0);
  x(0)=-1.0; x(N)=1.0;
  if (N==1) { return ; }

  VectorXd xint, w;
  JacobiGQ(alpha+1,beta+1,N-2,xint,w);

  // assemble result: sandwich eigenvalues between [-1,1]
  x(seq(1,N-1)) = xint;
  return ;
}
