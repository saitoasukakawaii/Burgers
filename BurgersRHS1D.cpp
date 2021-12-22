//
// Created by chen on 2021/11/30.
//
#include "burgers.h"
#include <cmath>

MatrixXd BurgersRHS1D(const double &time,
                  const MatrixXd &Dr, const MatrixXd &rx,
                  const MatrixXd &nx, const MatrixXd &LIFT,
                  const MatrixXd &Fscale, const MatrixXd &u) {
    MatrixXd du = MatrixXd::Constant(2, K, 0.);
    for (int i = 0; i < K - 1; i++) {
        du(1, i) = u(Nv - 1, i) - u(0, i + 1);
    }
    for (int i = 1; i < K; i++) {
        du(0, i) = u(0, i) - u(Nv - 1, i - 1);
    }
    double uin = -tanh((xmin + 0.5 - time) / (2 * epsilon)) + 1.0;
    du(0, 0) = 2.0 * (u(0, 0) - uin);
    double uout = -tanh((xmax + 0.5 - time) / (2 * epsilon)) + 1.0;
    du(1, K - 1) = 2.0 * (u(Nv - 1, K - 1) - uout);

    std::cout << "du:" << std::endl;
    std::cout << du << std::endl;
    MatrixXd dudx = rx.array()*(Dr * u).array();
    MatrixXd duB = 0.5*Fscale.array()*(nx.array()*du.array());
    MatrixXd q = sqrt(epsilon) * (dudx - LIFT * duB);
    std::cout << "q" << std::endl;
    std::cout << q << std::endl;
    MatrixXd dq = MatrixXd::Constant(2, K, 0.);
    for (int i = 0; i < K - 1; i++) {
        dq(1, i) = (q(Nv - 1, i) - q(0, i + 1)) / 2.;
    }
    for (int i = 1; i < K; i++) {
        dq(0, i) = (q(0, i) - q(Nv - 1, i - 1)) / 2.;
    }


    dq(0, 0) = 0.0;
    dq(1, K - 1) = 0.0;
    std::cout << "dq" << std::endl;
    std::cout << dq << std::endl;
    // Evaluate nonlinear flux
    MatrixXd du2 = MatrixXd::Constant(2, K, 0.);
    for (int i = 0; i < K - 1; i++) {
        du2(1, i) = ( u(Nv - 1, i)*u(Nv - 1, i) - u(0, i + 1)*u(0, i + 1) )/ 2.0;
    }
    for (int i = 1; i < K; i++) {
        du2(0, i) = (u(0, i)*u(0, i) - u(Nv - 1, i - 1)*u(Nv - 1, i - 1))/ 2.0;
    }
    // impose boundary condition
    du2(0, 0) = u(0, 0)*u(0, 0) - uin*uin;
    du2(1, K - 1) = u(Nv - 1, K - 1)*u(Nv - 1, K - 1) - uout*uout;
    std::cout << "du2" << std::endl;
    std::cout << du2 << std::endl;
    // Compute flux
    double maxvel = -1e5;
    for(int i=0;i<Nv;++i){
        for(int j=0;j<K;++j){
            if( abs(u(i,j))>maxvel ) maxvel = abs(u(i,j));
        }
    }
    std::cout << "maxvel" << std::endl;
    std::cout << maxvel << std::endl;
    //  penalty scaling--See Chapter 7.2
    // tau = .25 * reshape(N * N./ max(2 * J(vmapP), 2 * J(vmapM)), Nfp * Nfaces, K);
    double tau = 0;

    // flux term Lax - Friedrichs flux

    MatrixXd flux = MatrixXd(nx.array()*(0.5* du2  - sqrt(epsilon) * dq).array()) - 0.5*maxvel * du - sqrt(epsilon) * tau*du;
    std::cout << "flux" << std::endl;
    std::cout << flux << std::endl;
    // local derivatives of field
    MatrixXd dfdr = Dr * (MatrixXd(u.array().square() / 2) - sqrt(epsilon) * q);
    std::cout << "dfdr" << std::endl;
    std::cout << dfdr << std::endl;
    // compute right hand sides of the semi - discrete PDE
    MatrixXd rhsu;
    rhsu = -(MatrixXd(rx.array()*dfdr.array()) - LIFT * MatrixXd(Fscale.array()*flux.array()));
    std::cout << "rhsu" << std::endl;
    std::cout << rhsu << std::endl;
    return rhsu;
}