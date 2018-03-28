//
// Created by Michael Heuer on 28.03.18.
//

#ifndef AMOLQCPP_RADIALGAUSSIAN_H
#define AMOLQCPP_RADIALGAUSSIAN_H

#include <boost/math/special_functions.hpp>

class RadialGaussian{
public:

    RadialGaussian(double rCenter = 0, double sigma = 1/2.)
            : rCenter_(rCenter),
              sigma_(sigma),
              alpha_(1./(2*sigma*sigma)),
              normalizationConstant_g2_r2_(normalizationConstant_g2_r2())
    {}

    // computes the normalization constant for the integral S g^2 r^2 dr
    double normalizationConstant_g2_r2() const{
        double w = 2*alpha_;
        double W0 = 2*alpha_*rCenter_;
        double integral_r2_g2_dr =
                1./(4.*pow(w, 5/2.))
                *exp(-w*rCenter_*rCenter_)
                *(2*sqrt(w)*W0
                  + sqrt(M_PI)
                    *exp(W0*W0/w)
                    *(w+2*W0*W0)
                    *boost::math::erfc<double>(-W0/sqrt(w))
                );
        return  1./sqrt(integral_r2_g2_dr);
    };

    // computes the normalization constant for the integral S g r^2 dr
    double normalizationConstant_g_r2() const{
        double w = alpha_;
        double W0 = alpha_*rCenter_;
        double integral_r2_g_dr =
                1./(4.*pow(w, 2.5))*exp(-w*rCenter_*rCenter_)*(
                        2*sqrt(w)*W0 +
                        sqrt(M_PI)*exp(W0*W0/w)*(w+2*W0*W0)*(
                                1 - boost::math::erf<double>(-W0/sqrt(w))
                        )
                );
        return  1./sqrt(integral_r2_g_dr);
    };

    double value(double r) const {
        return std::exp(-alpha_*r*r);
    };

    double normalizedValue(double r) const {
        return value(r)*normalizationConstant_g2_r2_;
    };

    double alpha() const { return alpha_; };
    double center() const { return rCenter_; };

private:
    double rCenter_,sigma_,alpha_;
    double normalizationConstant_g2_r2_;
};

#endif //AMOLQCPP_RADIALGAUSSIAN_H
