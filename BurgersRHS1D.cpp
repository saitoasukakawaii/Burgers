//
// Created by chen on 2021/11/30.
//
#include "burgers.h"
#include <cmath>

void BurgersRHS1D(const double &time,
                  const MatrixXd &Dr,const MatrixXd &rx,
                  const MatrixXd &nx,const MatrixXd &LIFT,
                  MatrixXd &u, MatrixXd &rhsu)
{
    MatrixXd du = MatrixXd::Constant(2,K,0.);
    for(int i=0;i<K-1;i++)
{
        du(1,i) = u(Nv-1, i)-u(0, i+1);
}
for(int i=1;i<K;i++)
{
du(0,i) = u(0, i)-u(Nv-1, i-1);
}
    double uin =-tanh((xmin+0.5-time)/(2*epsilon))+1.0;
    du(0,0) = 2.0*(u(0,0)-uin);
    double uout=-tanh((xmax+0.5-time)/(2*epsilon))+1.0;
    du(1, K-1) = 2.0*(u(Nv-1, K-1)-uout);

    MatrixXd q = sqrt(epsilon)*(rx.*(Dr*u) - LIFT*(Fscale.*(nx.*du/2.0)));
    MatrixXd dq = atrixXd::Constant(2,K,0.);
    dq(:) = (q(vmapM)-q(vmapP))/2.0;

    dq(0,0) = 0.0; dq(1, K-1) = 0.0;
    %% Evaluate nonlinear flux
    du2 = zeros(Nfp*Nfaces,K); du2(:) = (u(vmapM).^2-u(vmapP).^2)/2.0;

    %% impose boundary condition
    du2(mapI)=(u(vmapI).^2-uin.^2); du2(mapO)=(u(vmapO).^2-uout.^2);

    %% Compute flux
    maxvel = max(max(abs(u)));

    %% penalty scaling -- See Chapter 7.2
                                      %tau = .25*reshape(N*N./max(2*J(vmapP),2*J(vmapM)), Nfp*Nfaces, K);
    tau=0;

    %% flux term
    flux = nx.*(du2/2.0 - sqrt(epsilon)*dq) - maxvel/2.0.*du - sqrt(epsilon)*tau.*du; % Lax-Friedrichs flux

                                                                                                       %% local derivatives of field
    dfdr = Dr*(u.^2/2 - sqrt(epsilon)*q);

    %% compute right hand sides of the semi-discrete PDE
            rhsu = -(rx.*dfdr - LIFT*(Fscale.*flux));
}