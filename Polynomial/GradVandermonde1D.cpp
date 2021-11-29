//
// Created by chen on 2021/9/15.
//

#include "Jacobi1D.h"

void GradVandermonde1D(const int &N, const VectorXd & r,
                           MatrixXd& DVr)
{
    int Nc = r.size();
    DVr = MatrixXd::Constant(Nc, N+1, 0.);
    VectorXd P;
    for(int i=0;i<N+1;i++)
    {
        GradJacobiP(r, 0, 0, i, P);
        DVr.col(i) = P;
    }
    return ;
}

