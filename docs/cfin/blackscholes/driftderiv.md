#### Derivation of drift from forward term structure in Black-Scholes Model 

We assume the forward price being a driftless geometric brownian motion such as: 

$$
\frac{dF_t(T)}{F_t(T)} = \sigma_t\left(F_t(T), t\right) dW_t
$$ 

Thus the dynamic of the log-process $\ln F_t(T)$ is given by: 

$$
d\ln F_t(T) = -\frac{1}{2}\sigma_t\left(F_t(T)\right)^2 dt + \sigma_t\left(F_t(T), t\right)dW_t
$$

The relationship between the log-spot price and the log-forward price is given by: 

$$\ln S_t = \lim_{T \rightarrow t} \ln F_t(T)$$

Thus: 

$$
d\ln S_t = d\ln F_t(T)\mid_{T=t} + \frac{\partial \ln F_t(T)}{\partial T}\mid_{T=t} dt = \left[\frac{\partial \ln F_t(T)}{\partial T}\mid_{T=t} - \frac{1}{2} \sigma_t(S_t, t)^2\right] dt + \sigma_t(S_t, t) dW_t
$$

Under the risk-neutral measure $\mathbb{Q}$, the forward price with maturity $T$ is given by:

$$
F_t(T) = \mathbb{E}^\mathbb{Q}_t\left[S_T\right]
$$

The solution of the log-spot at time $T$ is given by: 

$$
\ln S_T = \ln S_t + \int_t^T \mu_s - \frac{1}{2}\sigma^2_s(S_s, s)ds + \int_t^T \sigma_s(S_s, s) dW_s
$$

Thus: 

$$
S_T = S_t \exp \left[ \int_t^T \mu_s - \frac{1}{2}\sigma^2_s(S_s,s)ds + \int_t^T \sigma_s(S_s,s) dW_s \right]

$$

Finally we have: 

$$
F_t(T) = \mathbb{E}^\mathbb{Q}_t\left[S_T\right] = S_t\exp \left[ \int_t^T \mu_s ds\right] \text{ where } \mu_t = \frac{\partial \ln F_t(T)}{\partial T}\mid_{T=t}
$$