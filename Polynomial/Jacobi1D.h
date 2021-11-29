//
// Created by chen on 2021/9/18.
//

#ifndef DG_1D_BURGERS_JACOBI1D_H
#define DG_1D_BURGERS_JACOBI1D_H
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using namespace Eigen;

inline double SQ(double x){return x*x;}
class Jacobi1D
{
public:
    const int N;
    const int Np;
    VectorXd r;
    MatrixXd V, invV, Vr, Dr, invM;
    Jacobi1D(const int &_N);
    Jacobi1D(const Jacobi1D &Jac);
    ~Jacobi1D() = default;

};
void JacobiP(const VectorXd& x, const double &alpha, const double &beta, const int &N,
             VectorXd& P);
void Vandermonde1D(const int &N, const VectorXd& r,
                   MatrixXd& V1D);
void GradJacobiP(const VectorXd&,const double&,const double&,const int&,
                 VectorXd&);
void GradVandermonde1D(const int &N, const VectorXd & r,
                       MatrixXd&);
void JacobiGQ(const double&,const double&  ,const int&, // [in]
              VectorXd&   x,  VectorXd&   w);         // [out]
void JacobiGL(const double &alpha, const double &beta, const int &N,
              VectorXd& x);

#endif //DG_1D_BURGERS_JACOBI1D_H
