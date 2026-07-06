### Ordinary least squares 
The OLS model is a linear model given by: 

$$y = X\beta + \epsilon$$

where $y$ is a vector $(n\times1)$, $X$ is matrix $(n\times k)$, $\beta$ is a vector $(k\times 1)$ and $\epsilon$ is a vector $(n\times1)$. To fit the model we minimizes the sum of squared errors such that: 

$$\min_{\beta} (y - X\beta)^T(y - X\beta)$$

Thus setting to zero the first derivative with respect to $\beta$ gives: 

$$-2X\beta(y-X\beta) = 0 \Longleftrightarrow \hat{\beta} = (X^TX)^{-1}X^Ty$$

Plugging the true model within $\hat{\beta}$ gives: 

$$\hat{\beta} = \beta + (X^TX)^{-1}X^T\epsilon$$

This solution is only valid if the matrix $X^TX$ is invertible. If it's the case, we make the assumption that:

$$\epsilon \sim N(0, \sigma^2I)$$

where $\sigma^2$ is the true variance of the residuals and $I$ is the correlation matrix represented as the identity matrix. This implies: 
- Errors have a constant variance and no heteroskedasticity exist 
- Errors are i.i.d. which means that autocorrelation is null
From the assumption of normality, it follows that: 

$$E[\hat{\beta} \mid X] = \beta$$

$$Var[\hat{\beta} \mid X] = (X^TX)^{-1}X^T(\sigma^2I)X(X^TX)^{-1} = \sigma^2(X^TX)^{-1}$$

Which means that: 

$$\hat{\beta} \sim N(\beta, \sigma^2(X^TX)^{-1})$$

However, the true error variance is unknown and therefore estimated from observed residuals $\hat{\epsilon}$. Note first that the predicted values of $y$ is given by: 

$$
\hat{y} = X\hat{\beta} = X(X^TX)^{-1}X^Ty
$$ 

where $P = X(X^TX)^{-1}X^T$ is called the projection matrix and is a squared matrix of size $n$. First note some important properties about the projection matrix: 
- $PX = X$
- $P^T = P$ since $(ABC)^T = C^TB^TA^T$ and $X^TX$ is by definition symmetric
- $P^2 = P$ since $(X(X^TX)^{-1}X^T)(X(X^TX)^{-1}X^T) = X(X^TX)^{-1}X^T$
- $tr(P) = k$ since $tr(AB)=tr(BA)$ and $(X^TX)^{-1}X^TX$ is an identity matrix of size $k$
Thus let's define the observed residuals as $\hat{\epsilon} = y - \hat{y} = y - Py = (I-P)y$, since $(I-P)X = X-PX = 0$, leading to $\hat{\epsilon} = (I-P)\epsilon$. Thus: 

$$\hat{\epsilon}^T\hat{\epsilon} = \epsilon^TM\epsilon$$

given the previously described properties of matrix $P$. We note that we can diagonalized the matrix $M=(I-P)$ as: 

$$M = Q \Lambda Q^T$$

where $rank(\Lambda) = n-k$ given the properties of $P$. Since $\epsilon \sim N(0,\sigma^2I)$ then $z = \epsilon/\sigma \sim N(0,I)$ and $u = Q^Tz \sim N(0,I)$. Thus: 

$$\frac{\hat{\epsilon}^T\hat{\epsilon}}{\sigma^2} = u^T\Lambda u = \sum^{n-k}_{i=1}u_i \sim \chi_{n-k}$$

Alternatively, note that since $\hat{\sigma}^2 = \frac{\hat{\epsilon}^T\hat{\epsilon}}{n-k}$ then $\frac{(n-k)\hat{\sigma}^2}{\sigma^2} \sim \chi_{n-k}$ and it follows that: 

$$\frac{\hat{\beta_j} - \beta_j}{\sigma\sqrt{(X^TX)_{jj}^{-1}}} \times \frac{1}{\sqrt{\frac{(n-k)\hat{\sigma}^2}{\sigma^2}/(n-k)}} = \frac{\hat{\beta_j} - \beta_j}{\hat{\sigma}\sqrt{(X^TX)_{jj}^{-1}}} \sim t_{n-k}$$


#### C++ snippet example 

```cpp
int main() {

    Eigen::VectorXd y(15); 
    y << 52.21,53.12,54.48,55.84,57.20,58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46;

    Eigen::MatrixXd X(15,2);
    double heights[] = {1.47,1.50,1.52,1.55,1.57,1.6,1.63,1.65,1.68,1.70,1.73,1.75,1.78,1.80,1.83};

    for(int i=0;i<15;i++){
        X(i,0) = heights[i];
        X(i,1) = heights[i]*heights[i];
    }

    numcpp::reg::OLS ols(y, X, true);
    Eigen::VectorXd se = ols.standardErrors();
    Eigen::VectorXd ts = ols.tstats();
    Eigen::VectorXd pv = ols.pValues();

    std::cout << ols.betas(0) << std::endl;
    // 128.8128
    std::cout << se(0) << std::endl;
    // 16.3083
    std::cout << ts(0) << std::endl;
    // 7.8986
    std::cout << pv(0) << std::endl;
    // 0.0
    std::cout << ols.rSquared << std::endl;
    // 0.9989
    std::cout << ols.adjustedRSquared() << std::endl;
    // 0.9987

}
```