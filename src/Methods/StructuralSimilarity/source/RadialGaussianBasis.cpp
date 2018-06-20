//
// Created by Michael Heuer on 20.04.18.
//

#include "RadialGaussianBasis.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions/bessel.hpp>
#include <GaussKronrodCartesianIntegration.h>
#include <cmath>

//#include <mathimf.h>
//#include <mkl.h>


ModifiedSphericalBessel1stKind::ModifiedSphericalBessel1stKind(int degree) :
        _degree(degree) {
    _in.reserve(degree);
    _din.reserve(degree);
}

void ModifiedSphericalBessel1stKind::evaluate(double r, bool differentiate) {
    _in.clear();
    _din.clear();

    _in = ModifiedSphericalBessel1stKind::eval(_degree, r);

    if (differentiate) {
        _din.resize(_degree+1, 0.);
        if (r < RADZERO) {
            _din[1] = 1./3.;
            //_din.push_back(0.0);
            //_din.push_back(1./3.);
            //for (int n = 2; n <= _degree; ++n) {
            //    _din.push_back(0.);
            //}
        }
        else {
            _din[0] = _in[1];
            for (int n = 1; n <= _degree; ++n) {
                _din[n] = _in[n-1] - (n+1.)/r*_in[n];
            }
            //_din.push_back( _in[1] );
            //for (int n = 1; n <= _degree; ++n) {
            //    _din.push_back( _in[n-1] - (n+1.)/r*_in[n] );
            //}
        }
    }
    return;
}

std::vector<double> ModifiedSphericalBessel1stKind::eval(int degree, double r) {
    std::vector<double> il;
    if (r < RADZERO) {
        il.push_back(1.);
        il.push_back(0.);
    }
    else {
        il.push_back(sinh(r)/r);
        il.push_back(cosh(r)/r - sinh(r)/(r*r));
    }
    for (int l = 2; l <= degree; ++l) {
        if (r < RADZERO) {
            il.push_back(0.);
        }
        else {
            if (il[l-1] < SPHZERO) il.push_back(0.);
            il.push_back( il[l-2] - (2*(l-1)+1)/r*il[l-1] );
        }
    }
    return il;
}

RadialGaussianBasis::RadialGaussianBasis()
        : basis_(createBasis()),
          Sab_(Sab(ExpansionSettings::Radial::nmax)),
          radialTransform_(calculateRadialTransform(Sab_))
{}

std::vector<Gaussian> RadialGaussianBasis::createBasis() {
    const auto& nmax = ExpansionSettings::Radial::nmax;
    const auto& lmax = ExpansionSettings::Angular::lmax;
    const auto& sigmaAtom = ExpansionSettings::Radial::sigmaAtom;


    std::vector<Gaussian> basis;
    double basisFunctionCenter = 0;

    switch (ExpansionSettings::Radial::basisType){
        case ExpansionSettings::Radial::BasisType::equispaced : {
            for (int i = 0; i < nmax; ++i) {
                basisFunctionCenter = (ExpansionSettings::Cutoff::radius*double(i)) /double(nmax);
                basis.emplace_back(basisFunctionCenter, sigmaAtom);
            }
            return basis;
        }
        case ExpansionSettings::Radial::BasisType::adaptive :{
            basisFunctionCenter = 0;
            double sigmaStride = 1/2.;
            
            for (int i = 0; i < nmax; ++i) {
                double sigmaAdaptive = std::sqrt( 4./(2.*lmax + 1) * pow(basisFunctionCenter,2)+ pow(sigmaAtom,2) );
                basis.emplace_back(basisFunctionCenter, sigmaAdaptive);
                basisFunctionCenter += sigmaStride*sigmaAdaptive;
            }

            //TODO Attention: the cutoff radius is changed here
            // add lock function?
                    //TODO CHANGE THIS !!!!!!!!!!!!!!!!!!!!
            ExpansionSettings::Cutoff::radius = (*basis_.end()).center(); //TODO check this: is the last basis function centered at the cutoff radius?
            return basis;
        }
    }

}


double RadialGaussianBasis::operator()(double r, unsigned n) const{
    const auto& nmax  = ExpansionSettings::Radial::nmax;
    assert(n > 0 && "The radial basis function index must be positive");
    assert(n <= nmax && "The radial basis function index must be smaller than or equal to nmax");

    Eigen::VectorXd hvec(nmax);

    for (int i = 0; i < nmax; ++i) {
        hvec(i) = basis_[i].g2_r2_normalizedValue(r); //TODO store this?
    }
    return radialTransform_.col(n-1).dot(hvec);
};

Eigen::MatrixXd RadialGaussianBasis::Smatrix() const{
    return Sab_;
};

Eigen::MatrixXd RadialGaussianBasis::radialTransform() const{
    return radialTransform_;
};

Eigen::MatrixXd RadialGaussianBasis::Sab(unsigned nmax) const{
    Eigen::MatrixXd S(nmax,nmax);

    for (int i = 0; i < nmax; ++i) {
        for (int j = 0; j < nmax; ++j) { // skip iterations
            double a = basis_[i].alpha();
            double b = basis_[j].alpha();
            double rCenterA = basis_[i].center();
            double rCenterB = basis_[j].center();

            double w = a+b;
            double W0 = a*rCenterA + b*rCenterB;
            double s = 1./(4.*pow(w, 2.5));
            s *= exp(-a*rCenterA*rCenterA-b*rCenterB*rCenterB);
            s *= 2.0*sqrt(w)*W0
                 + sqrt(M_PI)*exp(std::pow(W0,2)/w)*(w+2*std::pow(W0,2))
                 * erfc(-W0/sqrt(w)); // TODO which one is faster (with MKL)
                 //*boost::math::erfc<double>(-W0/sqrt(w));
            s *= basis_[i].normalizationConstant_g2_r2()*basis_[j].normalizationConstant_g2_r2();
            S(i,j) = s;
        }
    }
    return S;
};

Eigen::MatrixXd RadialGaussianBasis::calculateRadialTransform(const Eigen::MatrixXd &Sab){
    Eigen::LLT<Eigen::MatrixXd> lltOfSab(Sab);
    auto L = lltOfSab.matrixL();
    return Eigen::Inverse<Eigen::MatrixXd>(L);
}

//TODO old integrand function => delete if other works
//double operator()(double rvar) const {
//    double exp_ik = exp(-beta_ik * (rvar - rho_ik) * (rvar - rho_ik));
//    return rvar * rvar * boost::math::sph_bessel<double>(l, 2 * ai * ri * rvar) * exp_ik;
//}

// Compute integrals S r^2 dr i_l(2*ai*ri*r) exp(-beta_ik*(r-rho_ik)^2) //TODO what is the difference between r and ri
std::vector<double> RadialGaussianBasis::calculateIntegrals(double ai, double ri, double rho_ik,double beta_ik) const {

    auto lmax = ExpansionSettings::Angular::lmax;
    auto n_steps = ExpansionSettings::Radial::integrationSteps;

    double sigma_ik = sqrt(0.5/beta_ik);
    double r_min = rho_ik - 4*sigma_ik;
    double r_max = rho_ik + 4*sigma_ik;
    if (r_min < 0.) {
        r_max -= r_min;
        r_min = 0.;
    }
    double delta_r_step = (r_max-r_min)/(n_steps-1);
    int n_sample = 2*n_steps+1;
    double delta_r_sample = 0.5*delta_r_step;

    ModifiedSphericalBessel1stKind mosbest(lmax);

    // Compute samples ...
    Eigen::MatrixXd integrand_l_at_r = Eigen::MatrixXd::Zero(lmax+1,n_sample);
    for (int s = 0; s < n_sample; ++s) {
        double r_sample = r_min - delta_r_sample + s*delta_r_sample;
        double exp_ik = exp(-beta_ik*(r_sample-rho_ik)*(r_sample-rho_ik));
        mosbest.evaluate(2*ai*ri*r_sample, false);// EQ 32 zweiter tiel
        //for (int l = 0; l != L_plus_1; ++l) {
        for (unsigned l = 0; l != lmax+1; ++l) {
            integrand_l_at_r(l,s) =
                    r_sample*r_sample*
                    mosbest._in[l]*
                    exp_ik;
        }
    }

    // ... integrate (à la Simpson)
    std::vector<double> ints(lmax+1);
    for (int s = 0; s < n_steps; ++s) {
        for (unsigned l = 0; l  != lmax+1; ++l ) {
            ints[l] += delta_r_step/6.*(
                    integrand_l_at_r(l, 2*s)+
                    4*integrand_l_at_r(l, 2*s+1)+
                    integrand_l_at_r(l, 2*s+2)
            );
        }
    }

    return ints;
}


double RadialGaussianBasis::calculateIntegral(double ai, double ri,unsigned l, double rho_ik,double beta_ik) const {


    // Sample coordinates along r-axis
    double sigma_ik = sqrt(0.5/beta_ik);
    double r_min = rho_ik - 4*sigma_ik;
    double r_max = rho_ik + 4*sigma_ik;
    if (r_min < 0.) {
        r_max -= r_min;
        r_min = 0.;
    }

    auto integrandFunction = [beta_ik,rho_ik,ai,ri,l](double r) -> double
    {
        double exp_ik = exp(-beta_ik * (r - rho_ik) * (r - rho_ik));
        return r * r * exp_ik
               *jn(l, 2 * ai * ri * r); // TODO which one is faster (with MKL)
               //* boost::math::sph_bessel<double>(l, 2 * ai * ri * r);
    };



    GaussKronrod::Integrator<double> integrator(ExpansionSettings::Radial::integrationSteps);
    GaussKronrod::Integrator<double>::QuadratureRule quadratureRule = Eigen::Integrator<double>::GaussKronrod15;

    double integral = integrator.quadratureAdaptive(integrandFunction, r_min, r_max,
            ExpansionSettings::Radial::desiredAbsoluteError, ExpansionSettings::Radial::desiredRelativeError, quadratureRule);
    return integral;
};

double RadialGaussianBasis::computeCoefficient(unsigned n, unsigned l, double centerToNeighborDistance, double neighborSigma) const {
    //double neighborSigma = sigma0_;

    if (neighborSigma < ExpansionSettings::Radial::radiusZero) {
        //TODO temporary, just to make it work. Lookup original implementation for actual treatment
        return 0;//!? basis_[n].g2_r2_normalizedValue(neighborPosition.norm() * radialTransform_(n-1,l));//TODO
    } else {



        //TODO Optimize
        double ai = 1/(2. * pow(neighborSigma,2));
        double ri = centerToNeighborDistance;
        SphericalGaussian gi_sph(Eigen::Vector3d::Zero(), neighborSigma); // <- position should not matter, as only normalization used here
        double norm_g_dV_sph_i = gi_sph.getNormalizationConstant();

        // Prefactor (r-independent)
            double basisFunctionAlpha = basis_[n-1].alpha(); //ak
            double basisFunctionCenter = basis_[n-1].center(); //rk
            double norm_r2_g2_dr_rad_k = basis_[n-1].normalizationConstant_g2_r2();
            double beta_ik = ai + basisFunctionAlpha;
            double rho_ik = basisFunctionAlpha * basisFunctionCenter / beta_ik;
            double prefac =
                    4 * M_PI *
                    norm_r2_g2_dr_rad_k * norm_g_dV_sph_i *
                    exp(-ai * ri * ri) *
                    exp(-basisFunctionAlpha * pow(basisFunctionCenter,2) * (1 - basisFunctionAlpha / beta_ik)); // eq 32 bzw. 33

        return prefac * calculateIntegral(ai, ri, l, rho_ik, beta_ik);
    }
}

Eigen::MatrixXd RadialGaussianBasis::computeCoefficients(double centerToNeighborDistance, double neighborSigma) const {
    //double neighborSigma = sigma0_;

    const auto lmax = ExpansionSettings::Angular::lmax;
    const auto nmax = ExpansionSettings::Radial::nmax;

    //TODO just resize?
    Eigen::MatrixXd radialCoeffsGnl = Eigen::MatrixXd::Zero(nmax,lmax+1);

    if (neighborSigma < ExpansionSettings::Radial::radiusZero) {

        for (unsigned n = 0; n < nmax; ++n) {
            double gn_at_r = basis_[n].value(centerToNeighborDistance);
            for (unsigned l = 0; l <= lmax; ++l) {
                radialCoeffsGnl(n, l) = gn_at_r;
            }
        }
        //return {};//!? basis_[n].g2_r2_normalizedValue(neighborPosition.norm() * radialTransform_(n-1,l));//TODO
    } else {



        //TODO Optimize
        double ai = 1 / (2. * pow(neighborSigma, 2));
        double ri = centerToNeighborDistance;
        SphericalGaussian gi_sph(Eigen::Vector3d::Zero(),
                                 neighborSigma); // <- position should not matter, as only normalization used here
        double norm_g_dV_sph_i = gi_sph.getNormalizationConstant();

        //iterate over all basis functions
        for (unsigned n = 0; n < nmax; ++n) {

            // Prefactor (r-independent)
            double basisFunctionAlpha = basis_[n].alpha(); //ak
            double basisFunctionCenter = basis_[n].center(); //rk
            double norm_r2_g2_dr_rad_k = basis_[n].normalizationConstant_g2_r2();
            double beta_ik = ai + basisFunctionAlpha;
            double rho_ik = basisFunctionAlpha * basisFunctionCenter / beta_ik;
            double prefac =
                    4 * M_PI *
                    norm_r2_g2_dr_rad_k * norm_g_dV_sph_i *
                    exp(-ai * ri * ri) *
                    exp(-basisFunctionAlpha * pow(basisFunctionCenter, 2) *
                        (1 - basisFunctionAlpha / beta_ik)); // eq 32 bzw. 33

            auto integrals = calculateIntegrals(ai, ri, rho_ik, beta_ik);

            for (unsigned l = 0; l <= lmax; ++l) {
                radialCoeffsGnl(n, l) = prefac * integrals[l];
            }
        }
    }
    //TODO VALIDATE! we forgot this so far
    //radialCoeffsGnl_ = radialTransform_.transpose() * radialCoeffsGnl_; //CAREFUL WITH TRANSPOSE
    return radialTransform_.transpose() * radialCoeffsGnl;
}
