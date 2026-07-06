### Brent 

The brent algorithm is defined in the wiki page: https://en.wikipedia.org/wiki/Brent%27s_method. It uses three optimization algortihms: Bisection, Secant and Inverse Quadratic Interpolation. More info from the latter are given here: https://medium.com/@swag1of3/inverse-quadratic-interpolation-c54a2f2e823b. Th idea of the latter is to solve a polynomial approximating the inverse of the target function which leads to the recursion formula for the next guess of the root. The algortihm is very fast to converge, however, if one of the three values $f(x_{n-2})$, $f(x_{n-1})$ or $f(x_n)$ coincide, the Inverse Quadratic Interpolation fails. A pseudo algortihm reads: 

1) Given $a,b$ and $c = \frac{a+b}{2}$, compute $f(a)$, $f(b)$ and $f(c)$
2) if $f(a) \ne f(c)$ and $f(b) \ne f(c)$, compute $s$ using inverse quadratic interpolation else compute $s$ using secant 
3) if some specific conditions are satisfied, replace $s$ using the bisection method 
4) given the values of $f(s)$, replace a or b by s and iterate 

One can also checks the wrok from Indranil Ghosh: https://indrag49.github.io/Numerical-Optimization/

#### C++ snippet example 




