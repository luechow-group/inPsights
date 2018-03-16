//
// Created by Michael Heuer on 15.03.18.
//

#ifndef AMOLQCPP_SPHERICALHARMONICSRADIALBASISEXPANDER_H
#define AMOLQCPP_SPHERICALHARMONICSRADIALBASISEXPANDER_H

#include <sh/spherical_harmonics.h>
#include <LebedevSphericalIntegration/SpatialFunction.h>
#include <SphericalIntegrator.h>
#include "RadialBasis.h"

#include <utility>
#include <vector>

class Coefficients{
public:
    Coefficients(std::vector<std::vector<std::vector<double >>> coefficients)
            : coefficients_(std::move(coefficients)){}
    double get(int n, int l, int m) const{
        return coefficients_[n-1][l][m+l];
    }
private:
    std::vector<std::vector<std::vector<double >>> coefficients_;
};

class SphericalHarmonicsRadialBasisExpander : public SpatialFunction{
    friend class SphericalIntegrator;

public:
    SphericalHarmonicsRadialBasisExpander(unsigned nmax, unsigned lmax , double rCutoff)
            : radialBasis_(nmax, rCutoff),
              sphericalIntegrator_(SphericalIntegratorSettings::expansion(lmax,nmax,rCutoff)),
              lmax_(lmax),nmax_(nmax),l_(0),n_(0),m_(0),rCutoff_(rCutoff),
              fPtr_(nullptr) {};


    Coefficients coefficients(SpatialFunction & f){
        std::vector<std::vector<std::vector<double >>> coefficients;
        for (int n = 1; n <= nmax_; ++n) {
            coefficients.emplace_back(std::initializer_list<std::vector<double >>{});

            for (int l = 0; l <= lmax_; ++l) {
                coefficients[n-1].emplace_back(std::initializer_list<double>{});

                for (int m = -l; m <= +l; ++m) {
                    //std::cout << (n-1) << l << m << " " << coefficient(f,n,l,m)<< std::endl;
                    coefficients[n-1][l].emplace_back(coefficient(f,n,l,m));
                }
            }
        }
        //std::cout << std::endl;
        return coefficients;
    };

    // f might be an AtomicNeighborhoodDensityFunction
    double coefficient(SpatialFunction & f, int n, int l, int m){
        set_n(n);
        set_l(l);
        set_m(m);

        setFunctionPtr(f);
        return sphericalIntegrator_.integrate(*this, rCutoff_);
    };

private:
    // This method is handed to the spherical integrator via *this, it overrides spatial function
    double value(const Eigen::Vector3d & rvec) const override {
        int idx = n_-1;
        auto l1 = l_;
        auto m2 = m_;
        auto rvecunit= rvec.normalized();
        return radialBasis_(rvec.norm(),idx)*sh::EvalSH(l_,m_,rvecunit)* fPtr_->value(rvec);
    };

    void setFunctionPtr(const SpatialFunction &f){
        fPtr_ = const_cast<SpatialFunction*>(&f);
    };

    void set_n(int n){
        assert(n > 0 && "n must be greater than zero.");
        assert(n <= nmax_ && "n must be smaller than or equal to nmax.");
        n_ = n;
    };
    void set_l(int l){
        assert(l >= 0 && "l must be greater than or equal to zero.");
        assert(l <= lmax_ && "l must be smaller than or equal to lmax.");
        l_ = l;
    };
    void set_m(int m){
        assert(m >= -lmax_ && "m must be greater than or equal to -lmax.");
        assert(m <= +lmax_ && "m must be smaller than or equal to lmax.");
        m_ = m;
    };



    RadialBasis radialBasis_;
    SphericalIntegrator sphericalIntegrator_;

    int lmax_, nmax_, l_, n_, m_;
    double rCutoff_;
    SpatialFunction* fPtr_; // needed to access specialized f_->value

};

#endif //AMOLQCPP_SPHERICALHARMONICSRADIALBASISEXPANDER_H
