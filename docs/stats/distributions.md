### distributions 

#### Uniform 

The Uniform distribution assigns equal probability density over an interval `[a,b]`. The sampling uses built-in from C++17. See the wiki page for details: https://en.wikipedia.org/wiki/Continuous_uniform_distribution

#### Normal 

The Normal (Gaussian) distribution is a continuous probability distribution defined by two parameters $\mu$ and $\sigma$. The sampling uses the inverse sampling method, using the inverse CDF method and the sampling method from the uniform ditribution class. Note the error function is given by the implementation from Jäckel (github.com/vollib/lets_be_rational). The inplementation of the inverse CDF is given by the Acklam's algorithm (https://web.archive.org/web/20151030215612/http://home.online.no/~pjacklam/notes/invnorm/). For more informations, see the wiki page: https://en.wikipedia.org/wiki/Normal_distribution

#### Gamma 

The Gamma distribution is a continuous probability distribution defined by two parameters $\theta$ and $k$. Sampling of the Gamma distribution is done using the Marsaglia's simple transformation-rejection method. See the wiki page for more info: https://en.wikipedia.org/wiki/Gamma_distribution#Random_variate_generation

#### Chi-Squared 

The chi-quared distribution with parameter $\nu$ is a special case of the Gamma distributions with $k = \nu/2$ and $\theta = 2$. 

#### T-Student

The T-Student is a continuous probability distribution with parameter $\nu$. See the wikki page for more info: https://fr.wikipedia.org/wiki/Loi_de_Student 


#### C++ snippet example

```cpp
int main() {

    numcpp::stats::ChiSquaredDistribution cs{2.0};
    numcpp::stats::GammaDistribution g{2.0,2.0};
    numcpp::stats::TStudentDistribution t{2.0};
    numcpp::stats::NormalDistribution n{2.0,3.0};


    std::cout << t.pdf(3.0) << std::endl;
    // 0.02741012223
    std::cout << g.pdf(7.0) << std::endl;
    // 0.05284542099
    std::cout << cs.pdf(7.0) << std::endl;
    // 0.01509869171
    std::cout << n.cdf(5.0) << std::endl;
    // 0.841344746

}
```