# Solving Burger's equation with spectral/hp element method

1D burger's equation with viscous term can be expreessed as:

 
$$\frac{\partial u}{\partial t}+\frac{\partial}{\partial x} \left(\frac{1}{2}u^2\right) = \nu \frac{\partial^2 u}{\partial x^2}$$

We use finite element method (FEM) to solve the equation. Selecting test function $\phi$ giving weighted residual integral formula:

$$(\frac{\partial u}{\partial t},\phi)+(\frac{\partial}{\partial x} \left(\frac{1}{2}u^2\right),\phi) = (\nu \frac{\partial^2 u}{\partial x^2},\phi)$$
