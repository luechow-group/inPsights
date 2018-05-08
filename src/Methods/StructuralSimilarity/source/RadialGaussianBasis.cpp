//
// Created by Michael Heuer on 20.04.18.
//

#include "RadialGaussianBasis.h"
#include "Environment.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <unsupported/Eigen/MatrixFunctions>


RadialGaussianBasis::RadialGaussianBasis(const ExpansionSettings &settings)
        : s_(settings),
          basis_(createBasis(s_)),
          Sab_(Sab(s_.radial.nmax)),
          radialTransform_(calculateRadialTransform(Sab_))
{}

std::vector<Gaussian> RadialGaussianBasis::createBasis(ExpansionSettings& settings) {

    std::vector<Gaussian> basis;
    double basisFunctionCenter = 0;
    
    switch (settings.radial.basisType){
        case RadialGaussianBasisType::equispaced : {
            for (int i = 0; i < settings.radial.nmax; ++i) {
                basisFunctionCenter = (settings.radial.cutoffRadius*i)/double(s_.radial.nmax);
                basis.emplace_back(basisFunctionCenter, settings.radial.sigmaAtom);
            }
        }
        case RadialGaussianBasisType::adaptive :{
            basisFunctionCenter = 0;
            double sigmaStride = 1/2.;
            
            for (int i = 0; i < settings.radial.nmax; ++i) {
                double sigmaAdaptive = std::sqrt( 4/(2.*settings.angular.lmax + 1) * pow(basisFunctionCenter,2)+ pow(settings.radial.sigmaAtom,2) );
                basis.emplace_back(basisFunctionCenter, sigmaAdaptive);
                basisFunctionCenter += sigmaStride*sigmaAdaptive;
            }

            // Attention: the cutoff radius is changed here 
            settings.radial.cutoffRadius = (*basis_.end()).center(); //TODO check this: is the last basis function centered at the cutoff radius?
        }
    }
    return basis;
}


double RadialGaussianBasis::operator()(double r, unsigned n) const{
    assert(n > 0 && "The radial basis function index must be positive");
    assert(n <= s_.radial.nmax && "The radial basis function index must be smaller than or equal to nmax");

    Eigen::VectorXd hvec(s_.radial.nmax);

    for (int i = 0; i < s_.radial.nmax; ++i) {
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
            double s;
            s = 1./(4.*pow(w, 2.5));
            s *= exp(-a*rCenterA*rCenterA-b*rCenterB*rCenterB);
            s *= 2*sqrt(w)*W0
                 + sqrt(M_PI)*exp(std::pow(W0,2)/w)*(w+2*std::pow(W0,2))
                   *boost::math::erfc<double>(-W0/sqrt(w));
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
double RadialGaussianBasis::calculateIntegral(double ai, double ri,unsigned l, double rho_ik,double beta_ik) const {
    int n_steps = 100;

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
        return r * r * boost::math::sph_bessel<double>(l, 2 * ai * ri * r) * exp_ik;
    };

    GaussKronrod::Integrator<double> integrator(n_steps);
    GaussKronrod::Integrator<double>::QuadratureRule quadratureRule = Eigen::Integrator<double>::GaussKronrod15;
    // Define the desired absolute and relative errors.
    double desAbsErr = 0;
    double desRelErr = Eigen::NumTraits<double>::epsilon() * 50;
    double integral = integrator.quadratureAdaptive(integrandFunction, r_min, r_max, desAbsErr, desRelErr, quadratureRule);

    return integral;
};

double RadialGaussianBasis::computeCoefficient(unsigned n, unsigned l, double centerToNeighborDistance, double neighborSigma) const {
    //double neighborSigma = sigma0_;

    if (neighborSigma < ZeroLimits::radiusZero) {
        //TODO temporary, just to make it work. Lookup original implementation for actual treatment
        return 0;//!? basis_[n].g2_r2_normalizedValue(neighborPosition.norm() * radialTransform_(n-1,l));//TODO
    } else {

        //TODO Optimize
        double ai = 1/(2. * pow(neighborSigma,2));
        double ri = centerToNeighborDistance;
        SphericalGaussian gi_sph(Eigen::Vector3d::Zero(), neighborSigma); // <- position should not matter, as only normalization used here
        double norm_g_dV_sph_i = gi_sph.getNormalizationConstant();

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
                    exp(-basisFunctionAlpha * pow(basisFunctionCenter,2) * (1 - basisFunctionAlpha / beta_ik)); // eq 32 bzw. 33

        return prefac * calculateIntegral(ai, ri, l, rho_ik, beta_ik);
    }
}
