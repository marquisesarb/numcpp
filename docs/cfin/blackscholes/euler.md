#### Euler discretization for the Blac-Scholes SDE

Given the black scholes SDE:

$$
dS_t/S_t = \mu_t dt + \sigma(S_t,t)dW_t \Longleftrightarrow dx_t = (\mu_t - \sigma^2(x_t,t)/2) dt + \sigma(x_t,t)dW_t
$$

where $x_t = \log S_t$. the solution of the SDE is given by:

$$
x_{t+\Delta t} = x_t + \int_t^{t+\Delta t} (\mu_s - \sigma^2(x_s,s)/2) ds + \int_t^{t+\Delta t} \sigma(x_s,s) dW_s
$$

since, under local volatility assumption, both quantities $\mu_t$ and $\sigma(x_t,t)$ are known at time $t$ we can write: 

$$
x_{t+\Delta t} = x_t + (\mu_t - \sigma^2(x_t,t)/2) \int_t^{t+\Delta t} ds + \sigma(x_t,t)\int_t^{t+\Delta t} dW_s
$$

which leads to: 

$$
x_{t+\Delta t} = x_t + (\mu_t - \sigma^2(x_t,t)/2)\Delta t + \sigma(x_t,t)\sqrt{\Delta T}Z
$$

with $ Z \sim N(0,1)$. Therefore: 

$$
S_{t+\Delta t} = S_t\exp^{(\mu_t - \sigma^2(x_t,t)/2)\Delta t + \sigma(x_t,t)\sqrt{\Delta T}Z}
$$