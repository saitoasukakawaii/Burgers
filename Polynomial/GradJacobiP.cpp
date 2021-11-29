// GradJacobiP.m
// function [dP] = GradJacobiP(z, alpha, beta, N);
// 2007/06/06
//---------------------------------------------------------
#include "Jacobi1D.h"

//---------------------------------------------------------
void GradJacobiP(const VectorXd& z,const double& alpha,const double& beta,const int& N,
                 VectorXd& dP)
//---------------------------------------------------------
{
  // function [dP] = GradJacobiP(z, alpha, beta, N);
  // Purpose: Evaluate the derivative of the orthonormal Jacobi
  //	   polynomial of type (alpha,beta)>-1, at points x
  //          for order N and returns dP[1:length(xp))]
  int Nc = z.size();
  dP = VectorXd::Constant(Nc, 0.);
  if (0 == N) {
    return ;
  } else {
      JacobiP(z,alpha+1,beta+1, N-1, dP);
    dP = sqrt(N*(N+alpha+beta+1))*dP;
  }
  return ;
}
