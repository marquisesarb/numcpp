#pragma once 
#include <cmath>
#include <numcpp/stats/distributions/uniform.hpp>

namespace numcpp::stats {

    inline double acklamStandardGaussianInverseCdf(double p) {

        const double split1 = 0.425;
        const double split2 = 5.0;
        const double const1 = 0.180625;
        const double const2 = 1.6;

        const double A0 = 3.3871328727963666080E0;
        const double A1 = 1.3314166789178437745E+2;
        const double A2 = 1.9715909503065514427E+3;
        const double A3 = 1.3731693765509461125E+4;
        const double A4 = 4.5921953931549871457E+4;
        const double A5 = 6.7265770927008700853E+4;
        const double A6 = 3.3430575583588128105E+4;
        const double A7 = 2.5090809287301226727E+3;
        const double B1 = 4.2313330701600911252E+1;
        const double B2 = 6.8718700749205790830E+2;
        const double B3 = 5.3941960214247511077E+3;
        const double B4 = 2.1213794301586595867E+4;
        const double B5 = 3.9307895800092710610E+4;
        const double B6 = 2.8729085735721942674E+4;
        const double B7 = 5.2264952788528545610E+3;

        const double C0 = 1.42343711074968357734E0;
        const double C1 = 4.63033784615654529590E0;
        const double C2 = 5.76949722146069140550E0;
        const double C3 = 3.64784832476320460504E0;
        const double C4 = 1.27045825245236838258E0;
        const double C5 = 2.41780725177450611770E-1;
        const double C6 = 2.27238449892691845833E-2;
        const double C7 = 7.74545014278341407640E-4;
        const double D1 = 2.05319162663775882187E0;
        const double D2 = 1.67638483018380384940E0;
        const double D3 = 6.89767334985100004550E-1;
        const double D4 = 1.48103976427480074590E-1;
        const double D5 = 1.51986665636164571966E-2;
        const double D6 = 5.47593808499534494600E-4;
        const double D7 = 1.05075007164441684324E-9;

        const double E0 = 6.65790464350110377720E0;
        const double E1 = 5.46378491116411436990E0;
        const double E2 = 1.78482653991729133580E0;
        const double E3 = 2.96560571828504891230E-1;
        const double E4 = 2.65321895265761230930E-2;
        const double E5 = 1.24266094738807843860E-3;
        const double E6 = 2.71155556874348757815E-5;
        const double E7 = 2.01033439929228813265E-7;
        const double F1 = 5.99832206555887937690E-1;
        const double F2 = 1.36929880922735805310E-1;
        const double F3 = 1.48753612908506148525E-2;
        const double F4 = 7.86869131145613259100E-4;
        const double F5 = 1.84631831751005468180E-5;
        const double F6 = 1.42151175831644588870E-7;
        const double F7 = 2.04426310338993978564E-15;

        if (p<=0) return log(p);
        if (p>=1) return log(1-p);

        const double q = p-0.5;
        if (fabs(q) <= split1)
        {
            const double r = const1 - q*q;
            return q * (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0) /
                (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + 1.0);
        }
        else
        {
            double r = q<0.0 ? p : 1.0-p;
            r = sqrt(-log(r));
            double ret;
            if (r < split2)
            {
                r = r - const2;
                ret = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0) /
                    (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + 1.0);
            }
            else
            {
                r = r - split2;
                ret = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0) /
                    (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + 1.0);
            }
            return q<0.0 ? -ret : ret;
        }
    }

    inline double dint(const double x) { return( (x>0) ? floor(x) : -floor(-x) ); }

    inline double calErf(double x, const int jint) {

        static const double a[5] = { 3.1611237438705656,113.864154151050156,377.485237685302021,3209.37758913846947,.185777706184603153 };
        static const double b[4] = { 23.6012909523441209,244.024637934444173,1282.61652607737228,2844.23683343917062 };
        static const double c__[9] = { .564188496988670089,8.88314979438837594,66.1191906371416295,298.635138197400131,881.95222124176909,1712.04761263407058,2051.07837782607147,1230.33935479799725,2.15311535474403846e-8 };
        static const double d__[8] = { 15.7449261107098347,117.693950891312499,537.181101862009858,1621.38957456669019,3290.79923573345963,4362.61909014324716,3439.36767414372164,1230.33935480374942 };
        static const double p[6] = { .305326634961232344,.360344899949804439,.125781726111229246,.0160837851487422766,6.58749161529837803e-4,.0163153871373020978 };
        static const double q[5] = { 2.56852019228982242,1.87295284992346047,.527905102951428412,.0605183413124413191,.00233520497626869185 };

        static const double zero = 0.;
        static const double half = .5;
        static const double one = 1.;
        static const double two = 2.;
        static const double four = 4.;
        static const double sqrpi = 0.56418958354775628695;
        static const double thresh = .46875;
        static const double sixten = 16.;

        double y, del, ysq, xden, xnum, result;

        static const double xinf = 1.79e308;
        static const double xneg = -26.628;
        static const double xsmall = 1.11e-16;
        static const double xbig = 26.543;
        static const double xhuge = 6.71e7;
        static const double xmax = 2.53e307;
        
        y = std::fabs(x);
        if (y <= thresh) {

            ysq = zero;

            if (y > xsmall) {
                ysq = y * y;
            }

            xnum = a[4] * ysq;

            xden = ysq;

            for (int i__ = 1; i__ <= 3; ++i__) {
                xnum = (xnum + a[i__ - 1]) * ysq;
                xden = (xden + b[i__ - 1]) * ysq;
            }
            result = x * (xnum + a[3]) / (xden + b[3]);
            if (jint != 0) {
                result = one - result;
            }
            if (jint == 2) {
                result = exp(ysq) * result;
            }
            goto L800;

        } else if (y <= four) {
            xnum = c__[8] * y;
            xden = y;
            for (int i__ = 1; i__ <= 7; ++i__) {
                xnum = (xnum + c__[i__ - 1]) * y;
                xden = (xden + d__[i__ - 1]) * y;
            }
            result = (xnum + c__[7]) / (xden + d__[7]);
            if (jint != 2) {
                double d__1 = y * sixten;
                ysq = dint(d__1) / sixten;
                del = (y - ysq) * (y + ysq);
                d__1 = exp(-ysq * ysq) * exp(-del);
                result = d__1 * result;
            }

        } else {
            result = zero;
            if (y >= xbig) {
                if (jint != 2 || y >= xmax) {
                    goto L300;
                }
                if (y >= xhuge) {
                    result = sqrpi / y;
                    goto L300;
                }
            }
            ysq = one / (y * y);
            xnum = p[5] * ysq;
            xden = ysq;
            for (int i__ = 1; i__ <= 4; ++i__) {
                xnum = (xnum + p[i__ - 1]) * ysq;
                xden = (xden + q[i__ - 1]) * ysq;
            }
            result = ysq * (xnum + p[4]) / (xden + q[4]);
            result = (sqrpi - result) / y;
            if (jint != 2) {
                double d__1 = y * sixten;
                ysq = dint(d__1) / sixten;
                del = (y - ysq) * (y + ysq);
                d__1 = exp(-ysq * ysq) * exp(-del);
                result = d__1 * result;
            }
        }
        L300:
        if (jint == 0) {
            result = (half - result) + half;
            if (x < zero) {
                result = -(result);
            }
        } else if (jint == 1) {
            if (x < zero) {
                result = two - result;
            }
        } else {
            if (x < zero) {
                if (x < xneg) {
                    result = xinf;
                } else {
                    double d__1 = x * sixten;
                    ysq = dint(d__1) / sixten;
                    del = (x - ysq) * (x + ysq);
                    y = exp(ysq * ysq) * exp(del);
                    result = y + y - result;
                }
            }
        }
        L800:
        return result;
    }

    inline double erf(double x) {return calErf(x, 0);}

    inline double erfc(double x) {return calErf(x, 1);}

    inline double erfcx(double x) {return calErf(x, 2);}

    struct NormalDistribution final : public ContinuousUnivariateProbabilityDistribution {

        static constexpr double ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;   
        static constexpr double ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887; 
        double mu,sigma;

        NormalDistribution(double mu_, double sigma_): mu(mu_), sigma(sigma_) {}
        NormalDistribution(): mu(0.0), sigma(1.0) {}

        double pdf(double x) const override {
            double z =(x -mu)/sigma;
            return ONE_OVER_SQRT_TWO_PI*std::exp(-.5*z*z)/sigma;
        }

        double cdf(double x) const override {

            return .5 * (1.0 + erf(ONE_OVER_SQRT_TWO*(x-mu)/sigma));
        }

        double invCdf(double p) const {return mu+sigma*acklamStandardGaussianInverseCdf(p);}

        std::complex<double> cf(double t) const override {

            std::complex<double> i(0.0, 1.0); 
            return std::exp(i * t * mu - 0.5 * sigma * sigma * t * t);
        }

        double sample(std::mt19937& gen) const override {
            return invCdf(UniformDistribution{}.sample(gen));
        }

        double mean() const {return mu;}
        double variance() const {return sigma*sigma;}
        double skewness() const {return 0.0;}
        double kurtosis() const {return 0.0;}



    };


    


}