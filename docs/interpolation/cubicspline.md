### Cubic Spline 

For a set of observed data $\{(x_0,y_0), ..,(x_n,y_n)\}$ with strictly increasing $x_0<...<x_n$, the cubic spline defined local cubic polynomial for each interval $[x_i, x_{i+1}]$: 
$$S_i(x) = a_i + b_ih_i+c_ih_i^2+d_ih_i^3$$
with $h_i = (x-x_i)$ and $i=0,1,..,n-1$. Since there are 4 coefficients and $n$ intervals, there are $4n$ unknowns to solve. By carefully specifying the conditions of the piecewise function $S(x)$ we can defined $4n$ equations: 
- Interpolation condition ($2n$ equations): $S_i(x_i) = y_i,S_i(x_{i+1}) = y_{i+1}$
- $C^1$ condition ($n-1$ equations): $S'_{i-1}(x_i) = S'_i(x_i)$
- $C^2$ condition ($n-1$ equations): $S''_{i-1}(x_i) = S''_i(x_i)$
- Natural boundaries conditions (2 equations): $S''(x_0) = 0, S''(x_n) = 0$

From the $C^2$ condition we can define the term $M_i$ being the second derivative of $S(x)$ at $x_i$ such that: 

$$M_i = S''_(i-1)(x_i) = S''_i(x_i) = 2c_i \Longleftrightarrow c_i = M_i/2$$

Furthermore, we also know $M_{i+1}$ at $i$ since:

$$M_{i+1} = 2c_i + 6d_ih_i \Longleftrightarrow d_i = \frac{M_{i+1}-M_i}{6h_i}$$

From the first condition we know that $S_i(x_i) = a_i \Longleftrightarrow y_i = a_i$. Furthermore we have: 

$$S_i(x_{i+1}) - S_i(x_i) = y_{i+1}-y_i = b_ih_i + c_ih^2_i + d_ih^3_i \Longleftrightarrow b_i = \frac{y_{i+1}-y_i}{h_i} - \frac{h_i}{6}(2M_i + M_{i+1})$$

Finally, from the $C^1$ conditions, we get: 

$$b_i = b_{i-1} + c_{i-1}h_{i-1} + 3d_{i-1}h_{i-1}^2 \Longleftrightarrow h_{i-1}M_{i-1} + 2(h_{i-1} + h_i)M_i + h_iM_{i+1} = 6 \left[\frac{y_{i+1}-y_i}{h_i} -\frac{y_i - y_{i-1}}{h_{i-1}}\right]$$

with $M_0 = 0$ and $M_n = 0$. This is a tri-diagonal system which can be efficiently solved using Thomas algorithm ($O(n)$ linear time complexity). Once the vector of $M$ is found, all parameters of the cubic spline can be implied. 
