// Copyright (C) 2018-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <unsupported/Eigen/MatrixFunctions>
#include <SpecialMathFunctions/ModifiedSphericalBessel1stKind.h>
#include <ErrorHandling.h>
#include <NaturalConstants.h>
#include <ErrorHandling.h>
#include "RadialBasis.h"
#include "SOAPSettings.h"

using namespace SOAP;

RadialBasis::RadialBasis()
        : basis_(createBasis()),
          Sab_(Sab(Radial::settings.nmax())),
          radialTransform_(calculateRadialTransform(Sab_))
{}

std::vector<Gaussian> RadialBasis::createBasis() {
    const auto nmax = Radial::settings.nmax();
    const auto lmax = Angular::settings.lmax();
    const auto sigmaAtom = Radial::settings.sigmaAtom();

    assert(nmax > 1 && "nmax must be greater than 1");

    std::vector<Gaussian> basis;
    double basisFunctionCenter = 0;

    switch (Radial::settings.basisType()){
        case Radial::BasisType::equispaced : {
            for (unsigned i = 0; i < nmax; ++i) {
                basisFunctionCenter = (Cutoff::settings.radius()*double(i)) /double(nmax-1);
                basis.emplace_back(basisFunctionCenter, sigmaAtom);
            }
            return basis;
        }
        case Radial::BasisType::adaptive :{
            //throw NotImplemented();

            basisFunctionCenter = 0;
            double sigmaStride = 1/2.;

            for (unsigned i = 0; i < nmax; ++i) {
                double sigmaAdaptive = std::sqrt( 4./(2.*lmax + 1) * std::pow(basisFunctionCenter,2)+ std::pow(sigmaAtom,2) );
                basis.emplace_back(basisFunctionCenter, sigmaAdaptive);
                basisFunctionCenter += sigmaStride*sigmaAdaptive;
            }

            //TODO Attention: the cutoff radius is changed here,
            // this can interfere with the global variable we use in ExpansionSettings::Radial::cutoff
            Cutoff::settings.radius = (*basis_.end()).center(); //TODO check this: is the last basis function centered at the cutoff radius?
            return basis;
        }
        default:
            throw NotImplemented();
    }

}


double RadialBasis::operator()(double r, unsigned n) const{
    const auto& nmax  = Radial::settings.nmax();
    assert(n > 0 && "The radial basis function index must be positive");
    assert(n <= nmax && "The radial basis function index must be smaller than or equal to nmax");

    Eigen::VectorXd hvec(nmax);

    for (unsigned i = 0; i < nmax; ++i) {
        hvec(i) = basis_[i].g2_r2_normalizedValue(r);
    }
    return radialTransform_.col(n-1).dot(hvec);
};

Eigen::MatrixXd RadialBasis::Smatrix() const{
    return Sab_;
};

Eigen::MatrixXd RadialBasis::radialTransform() const{
    return radialTransform_;
};

Eigen::MatrixXd RadialBasis::Sab(unsigned nmax) const{
    Eigen::MatrixXd S(nmax,nmax);

    for (unsigned i = 0; i < nmax; ++i) {
        for (unsigned j = 0; j < nmax; ++j) { // skip iterations
            double a = basis_[i].alpha();
            double b = basis_[j].alpha();
            double rCenterA = basis_[i].center();
            double rCenterB = basis_[j].center();

            double w = a+b;
            double W0 = a*rCenterA + b*rCenterB;
            double s = 1./(4.*std::pow(w, 2.5));
            s *= std::exp(-a*rCenterA*rCenterA-b*rCenterB*rCenterB);
            s *= 2.0*std::sqrt(w)*W0
                 + std::sqrt(Constant::pi)*std::exp(std::pow(W0,2)/w)*(w+2*std::pow(W0,2))
                 * std::erfc(-W0/std::sqrt(w)); // TODO which one is faster (with MKL)
                 //*boost::math::erfc<double>(-W0/sqrt(w));
            s *= basis_[i].normalizationConstant_g2_r2()*basis_[j].normalizationConstant_g2_r2();
            S(i,j) = s;
        }
    }
    return S;
};

Eigen::MatrixXd RadialBasis::calculateRadialTransform(const Eigen::MatrixXd &Sab){
    Eigen::LLT<Eigen::MatrixXd> lltOfSab(Sab);
    auto L = lltOfSab.matrixL();
    return Eigen::Inverse<Eigen::MatrixXd>(L);
}

/* Copied code from soapxx */
// Compute integrals S r^2 dr i_l(2*ai*ri*r) exp(-beta_ik*(r-rho_ik)^2) //TODO what is the difference between r and ri
std::vector<double> RadialBasis::calculateIntegrals(double ai, double ri, double rho_ik,double beta_ik) const {
    // Thresholds are defined to spare unnecessary computations.
    // Here, we use the same thresholds as https://github.com/capoe/soapxx/blob/master/src/soap/functions.hpp
    const auto comparisionEpsilon = General::settings.comparisonEpsilon();
    const double sphZeroThreshold = 1e-4;

    const auto lmax = Angular::settings.lmax();
    const auto n_steps = Radial::settings.integrationSteps();

    double sigma_ik = sqrt(0.5/beta_ik);
    double r_min = rho_ik - 4*sigma_ik;
    double r_max = rho_ik + 4*sigma_ik;
    if (r_min < 0.) {
        r_max -= r_min;
        r_min = 0.;
    }
    double delta_r_step = (r_max-r_min)/(n_steps-1);
    unsigned n_sample = 2*n_steps+1;
    double delta_r_sample = 0.5*delta_r_step;

    std::vector<double> modifiedSphericalBessel1stKindResults(lmax);


    // Compute samples ...
    Eigen::MatrixXd integrand_l_at_r = Eigen::MatrixXd::Zero(lmax+1,n_sample);
    for (unsigned s = 0; s < n_sample; ++s) {
        double r_sample = r_min - delta_r_sample + s*delta_r_sample;
        double exp_ik = std::exp(-beta_ik*(r_sample-rho_ik)*(r_sample-rho_ik));
        // EQ 32 second part
        modifiedSphericalBessel1stKindResults =
                ModifiedSphericalBessel1stKind::evaluateToMaxDegree(lmax, 2*ai*ri*r_sample,
                                                                    comparisionEpsilon, sphZeroThreshold);

        for (unsigned l = 0; l <= lmax; ++l)
            integrand_l_at_r(l,s) = std::pow(r_sample,2) * modifiedSphericalBessel1stKindResults[l] * exp_ik;
    }

    // ... integrate (à la Simpson)
    std::vector<double> ints(lmax+1);
    for (unsigned s = 0; s < n_steps; ++s) {
        for (unsigned l = 0; l  != lmax+1; ++l ) {
            ints[l] += delta_r_step/6.*(
                    integrand_l_at_r(l, 2*s)
                    + 4*integrand_l_at_r(l, 2*s+1)
                    + integrand_l_at_r(l, 2*s+2)
            );
        }
    }

    return ints;
}

/* Copied code from soapxx */
Eigen::MatrixXd RadialBasis::computeCoefficients(double centerToNeighborDistance, double neighborSigma) const {
    const auto lmax = Angular::settings.lmax();
    const auto nmax = Radial::settings.nmax();

    Eigen::MatrixXd radialCoeffsGnl = Eigen::MatrixXd::Zero(nmax,lmax+1);

    if (neighborSigma < Radial::settings.sigmaZeroThreshold()) {

        for (unsigned n = 0; n < nmax; ++n) {
            double gn_at_r = basis_[n].value(centerToNeighborDistance);
            for (unsigned l = 0; l <= lmax; ++l) {
                radialCoeffsGnl(n, l) = gn_at_r;
            }
        }
    } else {
        double ai = 1.0 / (2.0 * std::pow(neighborSigma, 2));
        double ri = centerToNeighborDistance;
        double norm_g_dV_sph_i = SphericalGaussian::calculateNormalizationConstant(neighborSigma);

        //iterate over all basis functions
        for (unsigned n = 0; n < nmax; ++n) {

            // Prefactor (r-independent)
            double basisFunctionAlpha = basis_[n].alpha(); //ak
            double basisFunctionCenter = basis_[n].center(); //rk
            double norm_r2_g2_dr_rad_k = basis_[n].normalizationConstant_g2_r2();
            double beta_ik = ai + basisFunctionAlpha;

            double prefac =
                    4.0 * Constant::pi *
                    norm_r2_g2_dr_rad_k * norm_g_dV_sph_i
                    * std::exp(-ai * ri * ri)
                    * std::exp(-basisFunctionAlpha
                               * std::pow(basisFunctionCenter, 2)
                               * (1.0 - basisFunctionAlpha / beta_ik)
                               ); // eq 32 bzw. 33

            double rho_ik = basisFunctionAlpha * basisFunctionCenter / beta_ik;
            auto integrals = calculateIntegrals(ai, ri, rho_ik, beta_ik);

            for (unsigned l = 0; l <= lmax; ++l) {
                radialCoeffsGnl(n, l) = prefac * integrals[l];
            }
        }
    }
    return radialTransform_.transpose() * radialCoeffsGnl;
}
