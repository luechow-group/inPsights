//
// Created by Michael Heuer on 15.03.18.
//

#ifndef AMOLQCPP_SPHERICALHARMONICSRADIALBASISEXPANDER_H
#define AMOLQCPP_SPHERICALHARMONICSRADIALBASISEXPANDER_H

#include <sh/spherical_harmonics.h>
#include <LebedevSphericalIntegration/SpatialFunction.h>
#include <SphericalIntegrator.h>
#include "RadialBasis.h"

#include <vector>

class SphericalHarmonicsRadialBasisExpander : public SpatialFunction{
    friend class SphericalIntegrator;

public:
    SphericalHarmonicsRadialBasisExpander(unsigned lmax, unsigned nmax, double rCutoff)
            : fPtr_(nullptr),
              radialBasis_(nmax, rCutoff),
              sphericalIntegrator_(SphericalIntegratorSettings::expansion(lmax,nmax,rCutoff)),
              lmax_(lmax),nmax_(nmax),l_(0),n_(0),m_(0),rCutoff_(rCutoff) {};


    std::vector<std::vector<std::vector<double >>> coefficients(SpatialFunction & f){
        std::vector<std::vector<std::vector<double >>> coefficients;
        coefficients = {};
        for (unsigned n = 1; n <= nmax_; ++n) {
            //coefficients[n-1].emplace_back({});
            for (unsigned l = 0; l < lmax_; ++l) {
                //coefficients[n-1][l].emplace_back({});
                for (int m = -lmax_; m < +lmax_; ++m) {
                    coefficients[n-1][l].emplace_back(coefficient(f,n,l,m));
                }
            }
        }
        return coefficients;
    };

    // f might be an AtomicNeighborhoodDensityFunction
    double coefficient(SpatialFunction & f, unsigned n, unsigned l, int m){
        set_n(n); set_l(l); set_m(m);

        setFunctionPtr(f);
        return sphericalIntegrator_.integrate(*this, rCutoff_);
    };

private:

    // This method is handed to the spherical integrator via *this, it overrides spatial function
    double value(const Eigen::Vector3d & rvec) const override {
        return radialBasis_(rvec.norm(),n_)*sh::EvalSH(l_,m_,rvec)* fPtr_->value(rvec);
    };

    void setFunctionPtr(const SpatialFunction &f){
        fPtr_ = const_cast<SpatialFunction*>(&f);
    };

    void set_n(unsigned n){
        assert(n > 0 && "n must be greater than zero.");
        assert(n <= nmax_ && "n must be smaller than or equal to nmax.");
        n_ = n;
    };
    void set_l(unsigned l){
        assert(l >= 0 && "l must be greater than or equal to zero.");
        assert(n <= lmax_ && "l must be smaller than or equal to lmax.");
        l_ = l;
    };
    void set_m(int m){
        assert(m >= -lmax_ && "m must be greater than or equal to -lmax.");
        assert(n <= +lmax_ && "m must be smaller than or equal to lmax.");
        m_ = m;
    };


    SpatialFunction* fPtr_; // needed to access specialized f_->value
    RadialBasis radialBasis_;
    SphericalIntegrator sphericalIntegrator_;


    unsigned lmax_, nmax_, l_, n_;
    int m_;
    double rCutoff_;

};

#endif //AMOLQCPP_SPHERICALHARMONICSRADIALBASISEXPANDER_H
