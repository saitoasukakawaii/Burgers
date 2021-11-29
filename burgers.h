//
// Created by chen on 2021/11/30.
//

#ifndef BURGERS_BURGERS_H
#define BURGERS_BURGERS_H
#include <Eigen/Dense>
using namespace Eigen;

const int N = 4;
const int Nv= N+1;
const int K = 8;
const double xmin = -1;
const double xmax =  1;
const double CFL  = 0.75;
const double epsilon   = 0.1;
const double FinalTime = 1.5;

namespace RK4{
    const unsigned int N_t = 5;
    double const rk4a[N_t] = { 0.0,
                               -567301805773.0/1357537059087.0,
                               -2404267990393.0/2016746695238.0,
                               -3550918686646.0/2091501179385.0,
                               -1275806237668.0/842570457699.0};

    double const rk4b[N_t] = {1432997174477.0/9575080441755.0,
                              5161836677717.0/13612068292357.0,
                              1720146321549.0/2090206949498.0,
                              3134564353537.0/4481467310338.0,
                              2277821191437.0/14882151754819.0};

    double const rk4c[N_t] = {0.0,
                              1432997174477.0/9575080441755.0,
                              2526269341429.0/6820363962896.0,
                              2006345519317.0/3224310063776.0,
                              2802321613138.0/2924317926251.0};
}

void MeshGen1D(const double& xmin, const double& xmax, const int& K, VectorXd& VX);
void GeometricFactors1D(const MatrixXd &x, const MatrixXd &Dr,
                        MatrixXd &rx, MatrixXd &J);
#endif //BURGERS_BURGERS_H
