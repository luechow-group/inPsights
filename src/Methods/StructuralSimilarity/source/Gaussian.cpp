//
// Created by Michael Heuer on 20.04.18.
//

#include<Gaussian.h>

IGaussian::IGaussian(double sigma = 1/2.)
        : sigma_(sigma), alpha_(1./(2*sigma*sigma)) {}


Gaussian::Gaussian(double rCenter = 0, double sigma = 1/2.)
        : IGaussian(sigma),
          rCenter_(rCenter),
          normalizationConstant_(normalizationConstant()),
          normalizationConstant_g2_r2_(normalizationConstant_g2_r2())
{}

// computes the normalization constant for the 1D volume integral
double Gaussian::normalizationConstant() const {
    return  sqrt(alpha_/M_PI);
};

double Gaussian::center() const { return rCenter_; };

// computes the normalization constant for the integral S g^2 r^2 dr
double Gaussian::normalizationConstant_g2_r2() const{
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
double Gaussian::normalizationConstant_g_r2() const{
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

double Gaussian::value(double r) const {
    return std::exp(-alpha_*(r-rCenter_)*(r-rCenter_));
};

double Gaussian::g2_r2_normalizedValue(double r) const {
    return value(r)*normalizationConstant_g2_r2_;
};


SphericalGaussian::SphericalGaussian(const Eigen::Vector3d& rCenter = Eigen::Vector3d::Zero(), double sigma = 1/2.)
        : IGaussian(sigma),
          rCenter_(rCenter),
          normalizationConstant_(normalizationConstant())
{};

// computes the normalization constant for the 3D volume integral
double SphericalGaussian::normalizationConstant() const {
    return pow(alpha_/M_PI, 3./2.);
};


double SphericalGaussian::value(const Eigen::Vector3d &r) const {
    return std::exp(-alpha_*(r-rCenter_).squaredNorm());
}

double SphericalGaussian::normalizedValue(const Eigen::Vector3d &r) const {
    return value(r) * normalizationConstant_;
}

const Eigen::Vector3d &SphericalGaussian::center() const { return rCenter_; };
