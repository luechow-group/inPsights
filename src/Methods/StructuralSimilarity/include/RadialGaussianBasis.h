//
// Created by Michael Heuer on 22.03.18.
//

#ifndef AMOLQCPP_RADIALGAUSSIANBASIS_H
#define AMOLQCPP_RADIALGAUSSIANBASIS_H

#include "Gaussian.h"
#include <Eigen/Core>
#include <boost/math/special_functions/bessel.hpp>
#include <GaussKronrodCartesianIntegration.h>

namespace ZeroLimits{
    const double radiusZero = 1e-10;
}

class RadialGaussianBasis{
public:
    // adaptive
    explicit RadialGaussianBasis(unsigned nmax, unsigned lmax, double sigma0 = 1 / 2.);

    // equispaced
    explicit RadialGaussianBasis(unsigned nmax, double rCut = 4.0, unsigned lmax = 4, double sigma = 1 / 2.); // sensible default?

    double operator()(double r, unsigned n) const;

    void computeCoefficients(const Eigen::Vector3d& rvec){ //TODO warum sigma als input benötigt?
        double particle_sigma = sigma0_;

        double r = rvec.norm();
        // Delta-type expansion =>
        // Second (l) dimension of <save_here> and <particle_sigma> ignored here
        if (particle_sigma < ZeroLimits::radiusZero) {

            for (unsigned n = 0; n < nmax_; ++n) {
                for (int l = 0; l < (lmax_+1)*(lmax_+1); ++l) {
                    //coefficients_(n, l) = basis_[n].value(r);
                    // use basis_[n].g2_r2_normalizedValue() instead?
                    coefficientsGnl_(n, l) = basis_[n].g2_r2_normalizedValue(r);
                }
            }
            //Gnl = ub::prod(_Tij, Gnl);
            // multiply with radial transform
            coefficientsGnl_ = radialTransform_.cwiseProduct(coefficientsGnl_);
        }
        else {
            // Particle properties
            double ai = 1. / (2 * particle_sigma * particle_sigma);
            double ri = r; // TODO Ist es richtig, rvec.norm() zu nehmen? Nomenklatur überdenken
            SphericalGaussian gi_sph(Eigen::Vector3d(0, 0, 0), particle_sigma); // <- position should not matter, as only normalization used here
            double norm_g_dV_sph_i = gi_sph.getNormalizationConstant();

            assert(basis_.size() == nmax_);

            for (unsigned k = 0; k < basis_.size(); ++k) {
                // Prefactor (r-independent)
                double ak = basis_[k].alpha();
                double rk = basis_[k].center();
                double norm_r2_g2_dr_rad_k = basis_[k].normalizationConstant_g2_r2();
                double beta_ik = ai + ak;
                double rho_ik = ak * rk / beta_ik;
                double prefac =
                        4 * M_PI *
                        norm_r2_g2_dr_rad_k * norm_g_dV_sph_i *
                        exp(-ai * ri * ri) *
                        exp(-ak * rk * rk * (1 - ak / beta_ik)); // eq 32 bzw. 33


                // Compute integrals S r^2 dr i_l(2*ai*ri*r) exp(-beta_ik*(r-rho_ik)^2)
                //std::vector<double> integrals;
                //integrals.resize(Gnl.size2(), 0.);
                //compute_integrals_il_expik_r2_dr(
                //        ai, ri, beta_ik, rho_ik, Gnl.size2(), _integration_steps,
                //        &integrals, NULL);
                for (unsigned l = 0; l < (lmax_+1)*(lmax_+1); ++l) {
                    coefficientsGnl_(k, l) = prefac * calculateIntegral(ai, ri, l, r, rho_ik, beta_ik);
                }
            }
        }
    };

    double getCoefficient(unsigned n, unsigned l){
        return coefficientsGnl_(n-1,l);
    }

    double sigma(){ return sigma0_; };

private:
    std::vector<Gaussian> createBasis(unsigned nmax, double rCut, double sigma);

    std::vector<Gaussian> createBasis(unsigned nmax, unsigned lmax, double sigma0);

    Eigen::MatrixXd Smatrix() const;

    Eigen::MatrixXd radialTransform() const;

    Eigen::MatrixXd Sab(unsigned nmax) const;

    Eigen::MatrixXd calculateRadialTransform(const Eigen::MatrixXd &Sab);


    //double operator()(double rvar) const {
    //    double exp_ik = exp(-beta_ik * (rvar - rho_ik) * (rvar - rho_ik));
    //    return rvar * rvar * boost::math::sph_bessel<double>(l, 2 * ai * ri * rvar) * exp_ik;
    //}

    // Compute integrals S r^2 dr i_l(2*ai*ri*r) exp(-beta_ik*(r-rho_ik)^2)
    double calculateIntegral(double ai, double ri,unsigned l, double r, double rho_ik,double beta_ik){
        int n_steps = 100;

        // Sample coordinates along r-axis
        double sigma_ik = sqrt(0.5/beta_ik);
        double r_min = rho_ik - 4*sigma_ik;
        double r_max = rho_ik + 4*sigma_ik;
        if (r_min < 0.) {
            r_max -= r_min;
            r_min = 0.;
        }

        /*
        double delta_r_step = (r_max-r_min)/(n_steps-1);
        int n_sample = 2*n_steps+1;
        double delta_r_sample = 0.5*delta_r_step;

        //ModifiedSphericalBessel1stKind mosbest(L_plus_1-1);

        // Compute samples ...
        ub::matrix<double> integrand_l_at_r = ub::zero_matrix<double>(L_plus_1, n_sample);
        for (int s = 0; s < n_sample; ++s) {
            double r_sample = r_min - delta_r_sample + s*delta_r_sample;
            double exp_ik = exp(-beta_ik*(r_sample-rho_ik)*(r_sample-rho_ik));
            //mosbest.evaluate(2*ai*ri*r_sample, gradients);// EQ 32 zweiter tiel
            for (int l = 0; l != L_plus_1; ++l) {
                integrand_l_at_r(l,s) =
                        r_sample*r_sample*
                        mosbest._in[l]*
                        boost::math::sph_bessel<double>(l,2*ai*ri*r_sample)*
                        exp_ik;
            }
        }
        // ... integrate (à la Simpson)
        std::vector<double> &ints = *integrals;
        for (int s = 0; s < n_steps; ++s) {
            //for (int l = 0; l != L_plus_1; ++l) {
                ints[l] += delta_r_step/6.*(
                        integrand_l_at_r(l, 2*s)+
                        4*integrand_l_at_r(l, 2*s+1)+
                        integrand_l_at_r(l, 2*s+2)
                );
            }
        }*/

        // integrand

        //auto f = [](double rvar,) {
        //    double exp_ik = exp(-beta_ik * (rvar - rho_ik) * (rvar - rho_ik));
        //    return rvar * rvar * boost::math::sph_bessel<double>(l, 2 * ai * ri * rvar) * exp_ik;
        //};


        auto integrandFunction = [beta_ik,rho_ik,ai,ri,l](double rvar) -> double
        {
            double exp_ik = exp(-beta_ik * (rvar - rho_ik) * (rvar - rho_ik));
            return rvar * rvar * boost::math::sph_bessel<double>(l, 2 * ai * ri * rvar) * exp_ik;
        };

        GaussKronrod::Integrator<double> integrator(n_steps);
        GaussKronrod::Integrator<double>::QuadratureRule quadratureRule = Eigen::Integrator<double>::GaussKronrod15;
        // Define the desired absolute and relative errors.
        double desAbsErr = 0;
        double desRelErr = Eigen::NumTraits<double>::epsilon() * 50;
        double integral = integrator.quadratureAdaptive(integrandFunction, r_min, r_max, desAbsErr, desRelErr, quadratureRule);

        return integral;
    };

    unsigned nmax_,lmax_;
    double sigma0_;
    std::vector<Gaussian> basis_;
    double rCut_;

    Eigen::MatrixXd Sab_, radialTransform_;
    Eigen::MatrixXd coefficientsGnl_;
};

#endif //AMOLQCPP_RADIALGAUSSIANBASIS_H
