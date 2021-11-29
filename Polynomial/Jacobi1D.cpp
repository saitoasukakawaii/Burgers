#include "Jacobi1D.h"


Jacobi1D::Jacobi1D(const int &_N): N(_N), Np(_N+1)
{
    JacobiGL(0,0,N,r);      // First get r
    Vandermonde1D(N,r,V);               // then get V
    GradVandermonde1D(N, r, Vr);        // then Vr
    invV = V.inverse();
    Dr = Vr*invV;
    invM = V*V.transpose();
}

Jacobi1D::Jacobi1D(const Jacobi1D &Jac): N(Jac.N), Np(Jac.Np)
{
    V = Jac.V;
    Vr = Jac.Vr;
    invV = Jac.invV;
    Dr = Jac.Dr;
    invM = Jac.invM;
}

