//
// Created by chen on 2021/11/30.
//

#include <iostream>
#include <fstream>
#include <sstream>
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

void output(const MatrixXd x, const MatrixXd u, const int nstep);

int main() {
    Jacobi1D jac(N);        // differential matrix and mass matrix
    VectorXd VX;            // vortexs
    MatrixXd x, J, rx;      // weight points of Gauss-quadrature
                            // and jacobian matrix

    // get the vortex coordinate of the 1D mesh
    double deltax = xmax - xmin;
    VX = VectorXd::Constant(K + 1, 0.);
    for (int i = 0; i < K + 1; i++) { VX(i) = deltax * i / K + xmin; }
    // get the physical coordinate of the weight points
    x = MatrixXd::Constant(N + 1, K, 0.);
    for (int i = 0; i < K; i++) {
        x(all, i) = 0.5 * ((1 - jac.r.array()) * VX(i) + (1 + jac.r.array()) * VX(i + 1));
    }
    // get the jacobian matrix
    J = jac.Dr * x;
    rx = 1. / J.array();
    // initial value of u
    MatrixXd u = -tanh((x.array() + 0.5) / (2 * epsilon)) + 1.0;
    double time = 0, timelocal;
    MatrixXd resu = MatrixXd::Constant(Nv, K, 0.);
    MatrixXd rhsu = MatrixXd::Constant(Nv, K, 0.);
    double umax = u.lpNorm<Infinity>();
    VectorXd dx = x(0, all) - x(1, all);
    double x_min = dx.array().abs().minCoeff();
    double dt = CFL * std::min(x_min / umax, x_min * x_min / std::sqrt(epsilon));
    int Nsteps = std::ceil(FinalTime / dt);
    dt = FinalTime / Nsteps;
    cout << "U max: " << umax << ", N step: " << Nsteps << ", delta t: " << dt << endl;
    // edge direction
    MatrixXd nx = MatrixXd::Constant(2, K, 0.);
    nx(0, all) = VectorXd::Constant(K, -1.0).transpose();   // negtive for left edge
    nx(1, all) = VectorXd::Constant(K, 1.0).transpose();    // positive for right edge
    MatrixXd LIFT = MatrixXd::Constant(Nv, 2, 0.);              // mass matrix of
    LIFT(all,0) = jac.invM(all,0);                       // left and right edge
    LIFT(all,1) = jac.invM(all,N);
    MatrixXd Fscale = MatrixXd::Constant(2, K, 0.);              // jacobian matrix of
    Fscale(0,all) = rx(0, all);                       // left and right edge
    Fscale(1,all) = rx(N, all);
////    output(x, u, 0);
    for (int tstep = 0; tstep< Nsteps; ++tstep) {
        for (int INTRK = 0; INTRK < RK4::N_t; ++INTRK) {
            timelocal = time + RK4::rk4c[INTRK] * dt;
            rhsu = BurgersRHS1D(timelocal, jac.Dr, rx, nx, LIFT, Fscale, u);
            resu = RK4::rk4a[INTRK] * resu + dt * rhsu;
            u = u + RK4::rk4b[INTRK] * resu;
        }
        // Increment time
        time = time + dt;
    }
//cout << rx(0,all)<<endl;
    return 0;
}

void output(const MatrixXd x, const MatrixXd u, const int nstep){
    std::ofstream ff;
    std::stringstream ss;
    ss << "U_" << std::setfill('0') << std::setw(8) << nstep << ".dat";
    ff.open(ss.str());
    for(int i=0;i<K;++i)
        for (int j = 0; j < Nv; ++j) {
            ff << x(j,i) << "\t" << u(j,i) << "\n";
        }
    ff.close();
}