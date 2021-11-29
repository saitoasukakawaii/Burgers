//
// Created by chen on 2021/11/30.
//
#include "burgers.h"
void GeometricFactors1D(const MatrixXd &x, const MatrixXd &Dr,
                        MatrixXd &rx, MatrixXd &J)
{

    MatrixXd xr = Dr*x;
    J = xr;
    int cols = J.cols();
    int rows = J.rows();
    rx.resize(rows, cols);
    rx= 1./J.array();
    return ;
}