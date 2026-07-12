#### Solving the Black-Scholes PDE

Under a log-normal local volatility (where Black-Scholes is a specific case of constant local volatility), $S_t$ satisfies the following geometric Brownian Motion SDE: 

$$
dS_t/S_t = \mu_t dt + \sigma(S_t,dt) dW_t
$$

where $\mu_t$ can be defined by the forward term structure $F_t$ of the spot process: 

$$
\mu_t = \frac{\partial \ln F_t(T)}{\partial T} \mid_{T=t}
$$

Furthermore, given a derivative asset $V_t$ written on $S_t$ and with boundary condition $V_T$, its dynamic can be written using Ito's lemma: 

$$
\partial_t V + \frac{\sigma(S,t)^2 S^2}{2} \partial^2_S V + \mu S \partial_S V - r V = 0
$$

By setting $x = \log(S/S_0)$ and using the chain-rule, one can obtain the simplified PDE: 

$$
\partial_t V + \frac{\sigma(x,t)^2}{2} \partial^2_x V + (\mu - \frac{1}{2}\sigma(x,t)^2) \partial_x V - r V = 0
$$

We directly make note here that the Crank Nicolson method for solving the PDE is just the average case of both the explicit and implicit methods. Thus we solely present both methods. 


####  <ins>__1) Explicit numerical method for solving the PDE__</ins>

The discretization is written as: 

$$
\frac{V^{j}_{n+1} - V_n^{j}}{\Delta T} + \frac{\sigma_l^2}{2} \frac{V^{j+1}_{n+1} + V^{j-1}_{n+1} - 2V^{j}_{n+1}}{(\Delta x_{n+1})^2} + (\mu - \frac{1}{2}\sigma_l^2) \frac{V^{j+1}_{n+1} - V^{j-1}_{n+1}}{2\Delta x_{n+1}} - rV^j_{n} = 0
$$

which leads to the following recursive scheme: 

$$
V^j_n = \frac{1}{1+r\Delta T} \left(p^{u}_{n+1} V^{j+1}_{n+1} + p^{m}_{n+1} V^{j}_{n+1} + p^{d}_{n+1} V^{j-1}_{n+1} \right)
$$

with: 

$$
p^{u} = \frac{\sigma_l^2 \Delta T}{2 (\Delta x)^2} + \frac{(\mu - \sigma_l^2/2)\Delta T}{2\Delta x}
$$
$$
p^{d} = \frac{\sigma_l^2 \Delta T}{2 (\Delta x)^2} - \frac{(\mu - \sigma_l^2/2)\Delta T}{2\Delta x}
$$
$$
p^{m} = 1 - \frac{\sigma_l^2 \Delta T}{(\Delta x)^2}
$$

where $\sigma_l = \sigma(x_j, t_n)$ is the local volatility point. it follows that $p^{u} + p^{d} + p^{m} = 1$. For stability, each quantity should behave as probability as this scheme is simply a trinominal tree. Note that the "probability" constraint requires that: 

$$ 
p^m \in (0,1) \Longleftrightarrow \frac{\sigma_l^2 \Delta T}{(\Delta x)^2} < 1 
$$

$$ 
p^u, p^d \in (0,1) \Longleftrightarrow \mid \mu - \sigma_l^2/2 \mid \Delta T / \Delta x < \frac{1}{2}
$$


####  <ins>__2) Implicit numerical method for solving the PDE__</ins>

The discretization is written as: 

$$
\frac{V^{j}_{n+1} - V_n^{j}}{\Delta T} + \frac{\sigma_l^2}{2} \frac{V^{j+1}_{n} + V^{j-1}_{n} - 2V^{j}_{n}}{(\Delta x)^2} + (\mu - \frac{1}{2}\sigma_l^2) \frac{V^{j+1}_{n} - V^{j-1}_{n}}{2\Delta x} - rV^j_{n} = 0
$$

When reorganized, this reads as: 

$$
V^{j}_{n+1} = \left(1 + r\Delta T + \frac{\sigma_l^2 \Delta T}{(\Delta x)^2}\right) V^j_n + \left(-\frac{\sigma_l^2 \Delta T}{2(\Delta x)^2} - \frac{(\mu - \sigma_l^2/2)\Delta T}{2\Delta x}\right) V^{j+1}_n + \left(-\frac{\sigma_l^2 \Delta T}{2(\Delta x)^2} + \frac{(\mu - \sigma_l^2/2)\Delta T}{2\Delta x}\right) V^{j-1}_n
$$

where $\sigma_l = \sigma(x_j, t_n)$ is the local volatility point.


####  <ins>__3) Payoff and their boundaries__</ins>


##### __Vanishing gamma boundary__

As the spot process converges such that $S \rightarrow \infty$, the payoff $V_t$ becomes linear: 

$$
\partial_x^2 V \rightarrow 0 \Longleftrightarrow \frac{V^M_{n} + V^{M-2}_{n} - 2V^{M-1}_{n}}{(\Delta x)^2} = 0 \Longleftrightarrow V^M_{n} = 2V^{M-1}_{n} - V^{M-2}_{n}
$$

##### __Absorbing state boundary__

The absorbing state $S_t = 0$ means that the payoff $V_t$ becomes fixed and display the following dynamic:

$$
\partial_t V - r V = 0 \Longleftrightarrow \frac{V^0_{n+1} - V^0_{n}}{\Delta T} - rV^0_{n} = 0 \Longleftrightarrow V^0_n = \frac{1}{1+ r\Delta T} V^0_{n+1}
$$

#####  __a) European and American vanilla payoffs__

Both european and american have the following terminal payoff $V_T$ written on $S_t$ with strike $K$: 

$$
V_T = \max(\theta(S_T - K), 0)
$$

where $\theta$ is equal to 1 for a call option and -1 for a put. Furthermore, since the american option has an early exercise feature, the holder will exercise if the price of the option today is less than the immediate payoff such that: 

$$
V_t = \max\left(\theta(S_t - K), \mathbb{E}_t\left[V_T\right]\right)
$$

