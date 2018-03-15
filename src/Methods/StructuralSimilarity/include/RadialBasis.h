//
// Created by Michael Heuer on 11.03.18.
//

#ifndef AMOLQCPP_RADIALBASIS_H
#define AMOLQCPP_RADIALBASIS_H

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

class RadialBasis{
public:
    RadialBasis(double rCutoff = 2, unsigned nmax = 4)
            : rCutoff_(rCutoff),
              W_(W(nmax)){
    };

    double NormalizationConstant(double rCutoff, double alpha){
        return std::sqrt( std::pow(rCutoff, 2*alpha+5) / (2*alpha+5) );
    }

    double phi(double r,double rCutoff, double alpha){
        return  std::pow( rCutoff-r, alpha+2) / NormalizationConstant(rCutoff,alpha);
    }

    Eigen::MatrixXd Sab(unsigned nmax){
        assert(nmax > 0 && "The number of radial basis functions must be greater than zero.");
        Eigen::MatrixXd Sab (nmax,nmax);

        for (unsigned i = 1; i <= Sab.rows(); ++i) {
            for (unsigned j = 1; j <= nmax; ++j) {
                Sab(i-1,j-1) = std::sqrt((5+2*i)*(5+2*j))/double(5+i+j);
            }
        }
        return Sab;
    }

    Eigen::MatrixXd W(unsigned nmax){
        assert(nmax > 0 && "The number of radial basis functions must be positive.");
        return Sab(nmax).inverse().sqrt();
    }

    double operator()(double r, unsigned i) {
        auto nMax = nmax();
        assert(i >= 0 && "The radial basis function index must be smaller than nmax");
        assert(i < nmax() && "The radial basis function index must be smaller than nmax");

        double sum = 0;
        for (int a = 0; a < nMax; ++a) {
            sum += W_(i,a) * phi(r, rCutoff_,a);
        }
        return sum;
    };

    unsigned long nmax(){
        assert(W_.rows() == W_.cols() && "The Matrix W must be square.");
        assert(W_.rows() > 0 && "Nmax must be positive.");
        return static_cast<unsigned long>(W_.rows());
    };

private:
    double rCutoff_;
    Eigen::MatrixXd W_;
};

#endif //AMOLQCPP_RADIALBASIS_H
