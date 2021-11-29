//
// Created by chen on 2021/11/30.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "Jacobi1D.h"
#include "burgers.h"
using std::endl;
using std::cout;
using std::string;
using namespace Eigen;



int main() {
    Jacobi1D jac(N);
    VectorXd VX;
    MatrixXd x, J, rx;
    MeshGen1D(xmin, xmax, K, VX);
    x = MatrixXd::Constant(N + 1, K, 0.);

    for (int i = 0; i < K; i++) {
        x(all, i) = 0.5 * ((1 - jac.r.array()) * VX(i) + (1 + jac.r.array()) * VX(i + 1));
    }

    GeometricFactors1D(x, jac.Dr, rx, J);
    MatrixXd u = -tanh((x.array() + 0.5) / (2 * epsilon)) + 1.0;
    double time = 0, timelocal;
    MatrixXd resu = MatrixXd::Constant(N + 1, K, 0.);
    MatrixXd rhsu = MatrixXd::Constant(N + 1, K, 0.);
    double umax = u.lpNorm<Infinity>();
    VectorXd dx = x(0, all) - x(1, all);
    double x_min = dx.array().abs().minCoeff();
    double dt = CFL * std::min(x_min / umax, x_min * x_min / std::sqrt(epsilon));
    int Nsteps = std::ceil(FinalTime / dt);
    dt = FinalTime / Nsteps;
    cout << umax << Nsteps << dt << endl;
    MatrixXd nx = MatrixXd::Constant(2, K, 0.);
    nx(0, all) = VectorXd::Constant(K, -1.0).transpose();
    nx(1, all) = VectorXd::Constant(K, 1.0).transpose();
    MatrixXd LIFT = MatrixXd::Constant(Nv, 2, 0.);
    LIFT(all,0) = jac.invM(all,0);
    LIFT(all,1) = jac.invM(all,N);


//    for (int tstep = 0; tstep< Nsteps; ++tstep){
//
//        for(int INTRK = 0; INTRK< RK4::N_t; ++INTRK) {
//            timelocal = time + RK4::rk4c[INTRK] * dt;
//            [rhsu] = BurgersRHS1D(u, epsilon, xL, xR, timelocal);
//            resu = RK4::rk4a[INTRK] * resu + dt * rhsu;
//            u = u + RK4::rk4b[INTRK] * resu;
//        }
//    // Increment time
//    time = time + dt
//    }

//    dt = CFL* min(xmin/umax,xmin^2/sqrt(epsilon));
//    Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;
//    for(int i=0;i)


    return 0;
}