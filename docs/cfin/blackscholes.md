### Black Scholes 

Implementation of the black scholes formula and its inversion for "implied volatility" with reference to Jäckel's work. For more info, see the following scientific papers: 
- [1] The pricing of options and corporate liabilities - Black and Scholes (1973) : https://www.jstor.org/stable/1831029
- [2] The pricing of commodity contract - Black (1976) : https://www.sciencedirect.com/science/article/abs/pii/0304405X76900246
- [3] Let's be rational - Jäckel (2015) : http://www.jaeckel.org/


The black scholes $S_t$ satisfies the following geometric Brownian Motion SDE: 

$$
dS_t/S_t = \mu_t dt + \sigma dW_t
$$

where $\mu_t$ is defined by the forward term structure $F_t$ of the spot process: 

$$
\mu_t = \frac{\partial \ln F_t(T)}{\partial T} \mid_{T=t}
$$

Furthermore, given a derivative asset $V_t$ written on $S_t$ and with boundary condition $V_T$, its dynamic can be written using Ito's lemma: 

$$
\partial_t V + \frac{\sigma^2 S^2}{2} \partial^2_S V + \mu S \partial_S V - r V = 0
$$

By setting $x = \log(S/S_0)$ and using the chain-rule, one can obtain the simplified PDE: 

$$
\partial_t V + \frac{\sigma^2}{2} \partial^2_x V + (\mu - \frac{1}{2}\sigma^2) \partial_x V - r V = 0
$$

#### 1) Explicit numerical method for solving the PDE

The discretization is written as: 

$$
\frac{V^{j}_{n+1} - V_n^{j}}{\Delta T} + \frac{\sigma^2}{2} \frac{V^{j+1}_{n+1} + V^{j-1}_{n+1} - 2V^{j}_{n+1}}{(\Delta x_{n+1})^2} + (\mu - \frac{1}{2}\sigma^2) \frac{V^{j+1}_{n+1} - V^{j-1}_{n+1}}{2\Delta x_{n+1}} - rV^j_{n} = 0
$$

which leads to the following recursive scheme: 

$$
V^j_n = \frac{1}{1+r\Delta T} \left(p^{u}_{n+1} V^{j+1}_{n+1} + p^{m}_{n+1} V^{j}_{n+1} + p^{d}_{n+1} V^{j-1}_{n+1} \right)
$$

with: 

$$
p^{u} = \frac{\sigma^2 \Delta T}{2 (\Delta x)^2} + \frac{(\mu - \sigma^2/2)\Delta T}{2\Delta x}
$$
$$
p^{d} = \frac{\sigma^2 \Delta T}{2 (\Delta x)^2} - \frac{(\mu - \sigma^2/2)\Delta T}{2\Delta x}
$$
$$
p^{m} = 1 - \frac{\sigma^2 \Delta T}{(\Delta x)^2}
$$

it follows that $p^{u} + p^{d} + p^{m} = 1$. For stability, each quantity should behave as probability as this scheme is simply a trinominal tree. Note that the "probability" constraint requires that: 

$$ 
p^m \in (0,1) \Longleftrightarrow \frac{\sigma^2 \Delta T}{(\Delta x)^2} < 1 
$$

$$ 
p^u, p^d \in (0,1) \Longleftrightarrow \mid \mu - \sigma^2/2 \mid \Delta T / \Delta x < \frac{1}{2}
$$

General boundaries condition can be applied which satisfies many different payoff. First note that the absorbing state $S_t = 0$ means that the payoff becomes fixed and display the following dynamic:

$$
\partial_t V - r V = 0 \Longleftrightarrow \frac{V^0_{n+1} - V^0_{n}}{\Delta T} - rV^0_{n} = 0 \Longleftrightarrow V^0_n = \frac{1}{1+ r\Delta T} V^0_{n+1}
$$

Similar analysis can be done when $S \rightarrow \infty$ as the payoff becomes linear: 

$$
\partial_x^2 V \rightarrow 0 \Longleftrightarrow \frac{V^M_{n} + V^{M-2}_{n} - 2V^{M-1}_{n}}{(\Delta x)^2} = 0 \Longleftrightarrow V^M_{n} = 2V^{M-1}_{n} - V^{M-2}_{n}
$$

Which leads to: 

$$
V^M_{n} =  \frac{1}{1+r\Delta T}\left(p^{u}_{n+1} V^M_{n+1} + (2p^{m}_{n+1} - p^u_{n+1}) V^{M-1}_{n+1} + (2p^{d}_{n+1} - p^m_{n+1}) V^{M-2}_{n+1} + p^{d}_{n+1}V^{M-3}_{n+1}\right)
$$



#### 2) Implicit numerical method for solving the PDE

The discretization is written as: 

$$
\frac{V^{j}_{n+1} - V_n^{j}}{\Delta T} + \frac{\sigma^2}{2} \frac{V^{j+1}_{n} + V^{j-1}_{n} - 2V^{j}_{n}}{(\Delta x)^2} + (\mu - \frac{1}{2}\sigma^2) \frac{V^{j+1}_{n} - V^{j-1}_{n}}{2\Delta x} - rV^j_{n} = 0
$$

When reorganized, this reads as: 

$$
V^{j}_{n+1} = \left(1 + r\Delta T + \frac{\sigma^2 \Delta T}{(\Delta x)^2}\right) V^j_n + \left(-\frac{\sigma^2 \Delta T}{2(\Delta x)^2} - \frac{(\mu - \sigma^2/2)\Delta T}{2\Delta x}\right) V^{j+1}_n + \left(-\frac{\sigma^2 \Delta T}{2(\Delta x)^2} + \frac{(\mu - \sigma^2/2)\Delta T}{2\Delta x}\right) V^{j-1}_n
$$

Similar to the explicit method, we have the following boundary conditions: 

$$
V^0_{n+1} = (1+ r\Delta T) V^0_n
$$

$$
V^M_{n+1} = 2V^{M-1}_{n+1} - V^{M-2}_{n+1}
$$

note that $V^M_{n+1}$ can be expressed in terms of $V^{M-1}_{n}$ and $V^{M-2}_{n}$