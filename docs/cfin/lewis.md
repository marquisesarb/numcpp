### Lewis Stochastic Volatiity framework
Generalization of the lewis framework for a given characteristic function using Gauss-Laguerre quadrature for numerical integration. For more info please check the foowing scientific paper from Lewis:

- [1] Option valuation under stochastic volatility - Lewis (2000) : https://econpapers.repec.org/bookchap/vsvvbooks/ovsv.htm
- [2] The volatility surface, a practitioner’s guide - Gatheral (2006) : https://link.springer.com/article/10.1007/s11408-007-0072-4

Pease note that the "normalized" european option price is implemented here folowing the formulation from P. Jäckel (Let's be rational). The derivation is given in the following section.

#### Derivation of lewis "normalized" european call price 

Given a future price process $\{F_t\}_{t\ge 0}$ and a call option $C$ written on the future with strike $K$, the covered call price at time $t$ is given by: 

$$
\Omega_t = F_t - C_t \Longleftrightarrow \Omega_T = \min(F_T,K)
$$ 

The normalized payoff $\omega$ is given by: 

$$\omega_t = \frac{\Omega_t}{\sqrt{F_tK}} = e^{x/2} - c_t = e^{x/2}\mathbb{E}_t[\min(e^{f_T}, e^{-x})]$$

with $x = \log(F_t/K)$ and $f_T = \log(F_T/F_t)$. $c_t$ denotes the normalized call price. The factorization from the min fonction can be done since $e^{x/2} > 0$. Before going further, we will recall the fourier transform G of a function $g$: 
$$
G(u) = \int_{-\infty}^{\infty} e^{iut} g(t) dt
$$
Also, note that the characteristic function of a random variable X defined as:
$$
\phi_X(u) = E[e^{iuX}] = \int_{\mathbb{R}} e^{iux} g_X(x)dx
$$
Thus, the characteristic function $\phi_X(t)$ is the fourier transform of the probability density function $g_X(x)$. Using this resut, one can find the fourier transofrm of the normalized covered payoff $\tilde{\omega}$ as:
$$
\tilde{\omega}_t(u) = \int_{\mathbb{R}} e^{iuk} e^{-k/2}\mathbb{E}_t[\min(e^{f_T},e^{k})] dk = e^{-k/2}\mathbb{E}_t\left[\int_{-\infty}^{f_T}e^{k(1+iu)}dk  + e^{f_T} \int_{f_t}^{v}e^{iuk}dk \right]
$$
with $k=-x$. Let us provide a sketch of proof for the converge of the complex exponential of the form $e^{iuk}$. Given $u \in \mathbb{C}$ and $a,b \in \mathbb{R}$, we decompose the complex number $u$ as $u = a + ib$. Euler provides an important identification for a complex exponential of the form $e^{ix}$ with $x \in \mathbb{R}$. The formula reads: 
$$
e^{ix} = \cos x + i \sin x
$$
This means the term $e^{ix}$ oscillates infinitely and does not have asymptotic limits on the real axis. Given the terms in $\tilde{\omega}(u,t)$, they can be decomposed as:

$$
e^{k(1+iu)} = e^{iak}e^{k(1-b)}
$$

$$
e^{iuk} = e^{iak}e^{-bk}
$$

Since $e^{iak}$ does not converge, the value of $b$ fully determines the global convergence of each exponential. Indeed, only when $b < 1$, one can write: 
$$
\lim_{k \rightarrow -\infty} e^{iak}e^{k(1-b)} = 0
$$
and only when $b > 0$, one can write: 
$$
\lim_{k \rightarrow \infty} e^{iak}e^{-bk}= 0
$$
Provided that the condition on $b = Im(u)$ holds, $\tilde{\omega}_t(u)$ can be simplified to: 

$$
\tilde{\omega}_t(u) = \frac{e^{x/2}}{u(u-i)}\phi_{f_T}(u-i)
$$

The inverse fourier transform allows to retrieve $g(t)$ from its fourier transform $G(u)$ using the following identity: 
$$
g(t) = \frac{1}{2\pi} \int_{-\infty}^{\infty} e^{-ist} G(s)ds
$$

Using this result, the covered call price $\omega_t$ can be expressed along the imaginary line $Im(u) = 1/2$. This "trick" (used by Gatheral in his book) insures that the fourier transform $\tilde{\omega}$ converges and allows for elegant algebraic simplification:
$$
\omega_t(x) = \frac{1}{2\pi} \int_{-\infty + i/2}^{\infty + i/2} e^{-iux} \tilde{\omega}_t(u)du =  \frac{1}{2\pi}\int_{-\infty}^{\infty} \frac{du}{u^2 + \frac{1}{4}}e^{iux}\phi_{f_T}(u-i/2) = \frac{1}{\pi}\int_{0}^{\infty} \frac{du}{u^2 + \frac{1}{4}}Re\left[e^{iux}\phi_{f_T}(u-i/2)\right]
$$


