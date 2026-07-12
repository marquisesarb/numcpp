#### Derivation of drift from forward term structure in Bachelier Model 

We assume the forward price being a driftless arithmetic brownian motion such as: 

$$
dF_t(T) = \sigma_t\left(F_t(T), t\right) dW_t
$$

Since the relationship between the spot price and the forward price is given by: 

$$ S_t = \lim_{T \rightarrow t} F_t(T)$$

Thus: 

$$
dS_t = dF_t(T)\mid_{T=t} + \frac{\partial F_t(T)}{\partial T}\mid_{T=t}dt = \frac{\partial F_t(T)}{\partial T}\mid_{T=t}dt + \sigma_t(S_t, t)dW_t
$$

With the solution given by: 

$$
S_T = S_t + \int_t^T \frac{\partial F_s(T)}{\partial T}\mid_{T=s}ds + \int_t^T \sigma_s(S_s, s)dW_s
$$

Finally, the forward price is given by: 

$$
F_t(T) = \mathbb{E}^\mathbb{Q}_t\left[S_T\right] = S_t + \int_t^T \mu_sds  \text{ where } \mu_t = \frac{\partial F_t(T)}{\partial T}\mid_{T=t} 
$$