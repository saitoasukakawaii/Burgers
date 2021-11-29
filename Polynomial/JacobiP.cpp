#include "Jacobi1D.h"

void JacobiP(const VectorXd& x, const double &alpha, const double &beta, const int &N,
                 VectorXd& P)
//---------------------------------------------------------
{
    // function [P] = JacobiP(x,alpha,beta,N)
    // Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
    //          (alpha+beta <> -1) at points x for order N and
    //          returns P[1:length(x)]
    // Note   : They are normalized to be orthonormal.
    // x: gauss point; N order of polynomial
    using boost::math::tgamma;
    double aold=0.0, anew=0.0, bnew=0.0, h1=0.0;
    double gamma0=0.0, gamma1=0.0;
    double ab=alpha+beta, ab1=alpha+beta+1.0, a1=alpha+1.0, b1=beta+1.0;

    int Nc = x.size();
    P = VectorXd::Constant(Nc, 0.);
    VectorXd  prow(Nc), x_bnew(Nc), Bnew(Nc), PROW(Nc), PP(Nc);
    MatrixXd PL(N+1, Nc);

    // Initial values P_0(x) and P_1(x)
    gamma0 = pow(2.0,ab1)/(ab1)*tgamma(a1)*tgamma(b1)/tgamma(ab1);
    double p0 = 1.0/sqrt(gamma0);
    P = VectorXd::Constant(Nc,p0);
    if (0==N) { return ; } else { PL.row(0) = P; }

    gamma1 = (a1)*(b1)/(ab+3.0)*gamma0;
    for (int i=0;i<Nc;i++){
        prow(i) = ((ab+2.0)*x(i)/2.0 + (alpha-beta)/2.0) / sqrt(gamma1);
    }
    P = prow;
    if (1==N) { return ; } else { PL.row(1) = prow; }

    // Repeat value in recurrence.
    aold = 2.0/(2.0+ab)*sqrt((a1)*(b1)/(ab+3.0));
    // Forward recurrence using the symmetry of the recurrence.
    for (int i=1; i<=(N-1); ++i) {
        h1 = 2.0*i+ab;
        anew = 2.0/(h1+2.0)*sqrt((i+1)*(i+ab1)*(i+a1)*(i+b1)/(h1+1.0)/(h1+3.0));
        bnew = - (SQ(alpha)-SQ(beta))/h1/(h1+2.0);
        Bnew = VectorXd::Constant(Nc, bnew);
        x_bnew = x-Bnew;
        for (int j=0;j<Nc;j++){prow(j) = x_bnew(j)*PL.row(i)(j);}
        PROW = -aold*PL.row(i-1);
        PP = PROW + prow;
        PL.row(i+1) = 1.0/anew*PP;
        aold = anew;
    }

    P = PL.row(N);
    return ;
}