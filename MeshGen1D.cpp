//
// Created by chen on 2021/11/30.
//
#include "burgers.h"
// K: the number of element, Nv the number of vortex point
void MeshGen1D(const double& xmin, const double& xmax, const int& K, VectorXd& VX)
{

    double dx = xmax-xmin;
    VX = VectorXd::Constant(K+1, 0.);
    for(int i=0;i<K+1;i++){VX(i)=dx*i/K+xmin;}
}
