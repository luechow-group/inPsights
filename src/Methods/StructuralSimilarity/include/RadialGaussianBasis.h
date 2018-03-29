//
// Created by Michael Heuer on 22.03.18.
//

#ifndef AMOLQCPP_RADIALGAUSSIANBASIS_H
#define AMOLQCPP_RADIALGAUSSIANBASIS_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/special_functions.hpp>
#include "RadialGaussian.h"

class RadialGaussianBasis{
public:
    explicit RadialGaussianBasis(int nmax, int lmax, double sigma0 = 1/2.)
            : nmax_(nmax),
              sigma0_(sigma0),
              basis_(createBasis(nmax,lmax,sigma0)),
              rCut_((*basis_.end()).center()),//TODO check this!
              Sab_(Sab(nmax)),
              radialTransform_(radialTransform(Sab_))
    {};

    explicit RadialGaussianBasis(int nmax, double rCut = 4.0, double sigma = 1/2.)
            : nmax_(nmax),
              sigma0_(sigma),
              basis_(createBasis(nmax,rCut,sigma)),
              rCut_(rCut),
              Sab_(Sab(nmax)),
              radialTransform_(radialTransform(Sab_))
    {};


    double operator()(double r, int n) const{
        assert(n > 0 && "The radial basis function index must be positive");
        assert(n <= nmax_ && "The radial basis function index must be smaller than or equal to nmax");

        Eigen::VectorXd hvec(nmax_);

        for (int i = 0; i < nmax_; ++i) {
            hvec(i) = basis_[i].normalizedValue(r); //TODO store this?
        }
        return radialTransform_.col(n-1).dot(hvec);
    };

    Eigen::MatrixXd Smatrix() const{
        return Sab_;
    }

    Eigen::MatrixXd radialTransform() const{
        return radialTransform_;
    }

private:
    std::vector<RadialGaussian> createBasis(int nmax, double rCut, double sigma) {

        std::vector<RadialGaussian> basis;
        for (int i = 0; i < nmax; ++i) {
            double rCenter = (rCut*i)/double(nmax);
            basis.emplace_back(rCenter,sigma);
        }
        return basis;
    };

    std::vector<RadialGaussian> createBasis(int nmax, int lmax, double sigma0) {
        double rCenter = 0;
        double sigmaStride = 1/2.;

        std::vector<RadialGaussian> basis;
        for (int i = 1; i <= nmax; ++i) {
            double sigma = std::sqrt( 4/(2.*lmax + 1) * rCenter*rCenter + sigma0*sigma0 );
            basis.emplace_back(rCenter,sigma);
            rCenter += sigmaStride*sigma;
        }
        return basis;
    };


    Eigen::MatrixXd Sab(int nmax) const{
        Eigen::MatrixXd S(nmax,nmax);

        for (int i = 0; i < basis_.size(); ++i) {
            for (int j = 0; j < basis_.size(); ++j) { // skip iterations
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


    Eigen::MatrixXd radialTransform(const Eigen::MatrixXd& Sab){
        Eigen::LLT<Eigen::MatrixXd> lltOfSab(Sab);
        auto L = lltOfSab.matrixL();
        return Eigen::Inverse<Eigen::MatrixXd>(L);
    }

    int nmax_;
    double sigma0_;
    std::vector<RadialGaussian> basis_;
    double rCut_;
    Eigen::MatrixXd Sab_, radialTransform_;
};

#endif //AMOLQCPP_RADIALGAUSSIANBASIS_H
