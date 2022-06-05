# Solving Burger's equation with discontinuous Galerkin method by spectral/hp element

1D burger's equation with viscous term can be expreessed as:

 
$$\frac{\partial u}{\partial t}+\frac{\partial}{\partial x} \left(\frac{1}{2}u^2\right) = \nu \frac{\partial^2 u}{\partial x^2}$$

We use finite element method (FEM) to solve the equation. Selecting test function $\phi$ giving weighted residual integral formula:

$$\left(\frac{\partial u}{\partial t},\phi\right)+\left(\frac{\partial}{\partial x} \left(\frac{1}{2}u^2\right),\phi\right) = \left(\nu \frac{\partial^2 u}{\partial x^2},\phi\right)$$

## Requirements
> Eigen-3.4

## How to use

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

Or just use awesome IDE, for example Clion to compile and run the program.
