//
// Created by Michael Heuer on 20.04.18.
//

#ifndef AMOLQCPP_SOAPEXPANSION_H
#define AMOLQCPP_SOAPEXPANSION_H

#include "RadialGaussianBasis.h"
#include "AngularBasis.h"
#include "CoefficientsVector.h"
#include <Particle.h>

//TODO find better name e.g. NeighborBasisExpander
template <typename Type>
class SoapExpansion{
public:
    explicit SoapExpansion(const ExpansionSettings& settings = ExpansionSettings::defaults())
            : s_(settings),
              radialGaussianBasis_(s_),
              angularBasis_(s_.angular)
    {}

    CoefficientsVector<Type> expandParticle(const Particle<Type> &center,
                                            const Particle<Type> &neighbor,
                                            double neighborSigma) const {
        double theta,phi;

        Eigen::Vector3d centerToNeighborVector = (neighbor.position()-center.position());
        double centerToNeighborDistance = centerToNeighborVector.norm();

        if(centerToNeighborDistance > 0.) {
            BoostSphericalHarmonics::ToSphericalCoords(centerToNeighborVector.normalized(), theta, phi);
            if (phi < 0.) phi += 2 * M_PI;
        } else { // center and neighbor positions are identical
            theta = 0.;
            phi = 0.;
        }

        CoefficientsVector<Type> coefficientsVector(s_);

        for (unsigned n = 1; n <= s_.radial.nmax; ++n) {
            for (unsigned l = 0; l <= s_.angular.lmax; ++l) {
                for (int m = -l; m <= l; ++m) {
                    auto coefficient = radialGaussianBasis_.computeCoefficient(n, l, centerToNeighborDistance, neighborSigma)
                                       * angularBasis_.computeCoefficient(l, m, theta, phi);

                    coefficientsVector.storeCoefficient(0, n, l, m, coefficient);
                }
            }
        }
        return coefficientsVector;
    }

    const ExpansionSettings& getSettings(){
        return s_;
    }

private:
    ExpansionSettings s_;
    RadialGaussianBasis radialGaussianBasis_;
    AngularBasis angularBasis_;
};

#endif //AMOLQCPP_SOAPEXPANSION_H
