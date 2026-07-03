### Newton-Raphson

See wiki : https://en.wikipedia.org/wiki/Newton%27s_method

The class NewtonRaphson has two constrcutors: 
- One constrcutor where both the function and its derivative must be separated 
- One constrcutor where a single function returns the value and its derivative through the std::pair C++17 interface 


The algortihm to modify the root's guess is given by: 

$$ x^{n+1} = x_n - \frac{f(x_n)}{f'(x_n)} $$

