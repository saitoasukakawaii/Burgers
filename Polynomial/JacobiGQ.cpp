
//---------------------------------------------------------
#include "Jacobi1D.h"
//---------------------------------------------------------
void JacobiGQ
(
  const double& alpha,  // [in]
  const double&  beta,   // [in]
  const int&     N,      // [in]
  VectorXd&   x,      // [out]
  VectorXd&   w      // [out]
)
//---------------------------------------------------------
{
  // function [x,w] = JacobiGQ(alpha,beta,N)
  // Purpose: Compute the N'th order Gauss quadrature points, x, 
  //          and weights, w, associated with the Jacobi 
  //          polynomial, of type (alpha,beta) > -1 ( <> -0.5).
    x = VectorXd::Constant(N+1, 0.);
    w = VectorXd::Constant(N+1, 0.);
  if (0==N) { x(0)=(alpha-beta)/(alpha+beta+2.0); w(0)=2.0; return; }

  double ab1 = alpha+beta+1.0, a1=alpha+1.0, b1=beta+1.0;
  MatrixXd Vr;
  MatrixXd J = MatrixXd::Constant(N+1,N+1,0.);


  //#######################################################
  // Note: this assembly differs from the Matlab script.
  //       - manual assembly of diagonals
  //       - sorting of LAPACK eigensystem
  //#######################################################
    VectorXd h1 = 2.0*VectorXd::LinSpaced(N+1,0,N)+VectorXd::Constant(N+1, alpha+beta);
    VectorXd d0(N+1);

  // main diagonal: diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1)
  double fac = -0.5*(SQ(alpha)-SQ(beta));
  const double eps = 1e-16;
  for (int i=0; i<=N; ++i) {
      if(i==0){
          // Handle division by zero
          if (alpha+beta < 10*eps) { d0(1)=0.0; }
          else { d0(i) = fac / (h1(i) + 2.0) / h1(i); }
      } else {
          d0(i) = fac / (h1(i) + 2.0) / h1(i);
      }
  }
    J = J+d0.asDiagonal().toDenseMatrix();
  // 1st upper diagonal: diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta) .* ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
  for (int i=1; i<=N; ++i) { J(i-1,i) = 2.0/(h1(i-1)+2.0)*sqrt(i*(i+alpha+beta)*(i+alpha)*(i+beta)/(h1(i-1)+1)/(h1(i-1)+3.0)); }

  // Form symmetric matrix from recurrence.
    MatrixXd J_ = J.transpose();
  J = J + J_;    // J = J + J';

  // Compute quadrature by eigenvalue solve
  // [Vr,D] = eig(J); x = diag(D);
    SelfAdjointEigenSolver<MatrixXd> eigensolver(J);
    if (eigensolver.info() != Success) abort();
    Vr = eigensolver.eigenvectors();
    x = eigensolver.eigenvalues();

//  if (sort) {
//
//    // Note: Matlab appears to sort results from eig()
//    // so that the eigenvalues are in ascending order.
//    // Here we sort the columns of eigenvector matrix
//    // with the same permutation required to sort the
//    // eigenvalues into ascending order. Target: make
//    // w=1st row of eigenvector matrix match Matlab's.
//
//    DVecSort sx=x;  IVec idx;
//    sx.makeIndex(idx);        // find sort permutation
//    sx.sortFromIndexVec(idx); // permute eigenvalues
//    Vr.sort_cols(idx);        // permute eigenvectors
//    x = sx;                   // copy sorted evals to x
//  }
  double coef = pow(2.0,ab1) / (ab1)*tgamma(a1)*tgamma(b1)/tgamma(ab1);
  for(int i=0;i<=N;i++){ w(i) = SQ(Vr(0,i))*coef;}
}
