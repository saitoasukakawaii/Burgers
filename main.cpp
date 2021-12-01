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



void output_result(const MatrixXd &x, const MatrixXd &u, const int &nstep);

int main() {
//    Jacobi1D jac(N);
//    cout.precision(15);
//    cout << std::fixed;
//    cout << jac.r;
//    VectorXd VX;
//    MatrixXd x, J, rx;
//    MeshGen1D(xmin, xmax, K, VX);
//    x = MatrixXd::Constant(N + 1, K, 0.);
//
//    for (int i = 0; i < K; i++) {
//        x(all, i) = 0.5 * ((1 - jac.r.array()) * VX(i) + (1 + jac.r.array()) * VX(i + 1));
//    }
//    std::cout.precision(15);
//    std::cout << std::fixed;
//    GeometricFactors1D(x, jac.Dr, rx, J);
//    cout << jac.Dr << rx << endl;
//    double u[Nv][K]={0.}, resu[Nv][K]={0.}, rhsu[Nv][K]={0}, LIFT[Nv][2]={0.};
//    double umax=0;
//    for(int i=0;i<K;++i)
//        for(int j=0;j<Nv;++j)
//        {
//            u[j][i] = -tanh((x(j,i) + 0.5) / (2 * epsilon)) + 1.0;
//            if(std::abs(u[j][i])>umax) umax = std::abs(u[j][i]);
//        }
//    double time = 0, timelocal;
//
//    VectorXd dx = x(0, all) - x(1, all);
//    double x_min = dx.array().abs().minCoeff();
//    double dt = CFL * std::min(x_min / umax, x_min * x_min / std::sqrt(epsilon));
//    int Nsteps = std::ceil(FinalTime / dt);
//    dt = FinalTime / Nsteps;
//
//
//    for(int j=0;j<Nv;++j)
//    {
//        LIFT[j][0] = jac.invM(j, 0);
//        LIFT[j][1] = jac.invM(j, N);
//    }
//
//    MatrixXd Fscale = MatrixXd::Constant(2, K, 0.);
//    Fscale(0,all) = rx(0, all);
//    Fscale(1,all) = rx(Nv-1, all);
//    output_result(x,u, 0);


//    cout << u(Nv - 1, 0)  << "\t" <<  u(0,  1) << endl;
//    for (int tstep = 0; tstep< 10; ++tstep){
//        std::cout << "\n\n";
//        std::cout << "----------------------------------------------------\n";
//        std::cout << "Time step: " << tstep << "\n" << std::endl;
//        for(int INTRK = 0; INTRK< RK4::N_t; ++INTRK) {
//            std::cout << "*******************************\n";
//            std::cout << "RK4 step: " << INTRK << "\n" << std::endl;
//            timelocal = time + RK4::rk4c[INTRK] * dt;
//            std::cout << "timelocal: " << timelocal << "\n" << std::endl;
//            BurgersRHS1D(timelocal,jac.Dr, rx, nx, LIFT, Fscale, u, rhsu);
//            resu = RK4::rk4a[INTRK] * resu + dt * rhsu;
//            std::cout << "resu: " << "\n" ;
//            std::cout << resu << "\n" ;
//            u = u + RK4::rk4b[INTRK] * resu;
//            std::cout << "u: "<< "\n" ;
//            std::cout << u << "\n" ;
//            std::cout << "*******************************\n";
//        }
//        // Increment time
//
//        time = time + dt;
//        output_result(x, u, tstep+1);
//        std::cout << "\n";
//        std::cout << "----------------------------------------------------\n";
//    }


    return 0;
}

void output_result(const MatrixXd &x, const MatrixXd &u, const int &nstep){
    std::stringstream tmp;
    tmp << "U_" << std::setfill('0') << std::setw(5) << nstep << ".dat";
    std::string Uname = tmp.str();
    std::ofstream ss(Uname);
    ss.precision(10);
    ss << std::fixed;
    for(int i=0;i<K;++i)
        for(int j=0;j<Nv;++j)
        {
            ss << x(j,i) << "\t" << u(j,i) << "\n";
        }
    ss.close();
}