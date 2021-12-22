# Solving Burger's equation with spectral/hp element method

1D burger's equation with viscous term can be expreessed as:

 
> <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mrow><mo>$</mo></mrow><mfrac><mrow><mi mathvariant="normal">∂</mi><mi>u</mi></mrow><mrow><mi mathvariant="normal">∂</mi><mi>t</mi></mrow></mfrac><mo>+</mo><mfrac><mi mathvariant="normal">∂</mi><mrow><mi mathvariant="normal">∂</mi><mi>x</mi></mrow></mfrac><mrow data-mjx-texclass="INNER"><mo data-mjx-texclass="OPEN">(</mo><mfrac><mn>1</mn><mn>2</mn></mfrac><msup><mi>u</mi><mn>2</mn></msup><mo data-mjx-texclass="CLOSE">)</mo></mrow><mo>=</mo><mi>ν</mi><mfrac><mrow><msup><mi mathvariant="normal">∂</mi><mn>2</mn></msup><mi>u</mi></mrow><mrow><mi mathvariant="normal">∂</mi><msup><mi>x</mi><mn>2</mn></msup></mrow></mfrac><mrow><mo>$</mo></mrow></math>

We use finite element method (FEM) to solve the equation. Selecting test function $\phi$ giving weighted residual integral formula:

> $(\frac{\partial u}{\partial t},\phi)+(\frac{\partial}{\partial x} \left(\frac{1}{2}u^2\right),\phi) = (\nu \frac{\partial^2 u}{\partial x^2},\phi)$
