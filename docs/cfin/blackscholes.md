### Black Scholes 

Implementation of the black scholes formula and its inversion for "implied volatility" with reference to Jäckel's work. For more info, see the following scientific papers: 
- [1] The pricing of options and corporate liabilities - Black and Scholes (1973) : https://www.jstor.org/stable/1831029
- [2] The pricing of commodity contract - Black (1976) : https://www.sciencedirect.com/science/article/abs/pii/0304405X76900246
- [3] Let's be rational - Jäckel (2015) : http://www.jaeckel.org/

#### Black scholes simulation 

Implementation of Euler discretization scheme for the Black Scholes SDE using a given forward term structure ($f : \mathbb{R_+} \rightarrow \mathbb{R}$). Given the existence of a forward curve $F_t(T)$ for the spot price $S_t$ such that $S_t = \lim_{T\rightarrow t} F_t(T)$, since the black scholes SDE satisfies a geometric brownian motion, the term structure function is given by: 
$$f(t) = \frac{\partial \ln F_t(T)}{\partial T} \mid_{T=t}$$
Constant, zero and maturity dependant drift can be applied within the simulation. 