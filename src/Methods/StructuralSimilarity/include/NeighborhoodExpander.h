//
// Created by Michael Heuer on 20.04.18.
//

#ifndef AMOLQCPP_NEIGHBORHOODEXPANDER_H
#define AMOLQCPP_NEIGHBORHOODEXPANDER_H

#include "RadialGaussianBasis.h"
#include "AngularBasis.h"
#include "NeighborhoodExpansionCoefficientsVector.h"
#include <Particle.h>

template <typename Type>
class NeighborhoodExpander{
public:
    explicit NeighborhoodExpander(unsigned numberOfParticles, const ExpansionSettings& settings = ExpansionSettings::defaults())
            : numberOfParticles_(numberOfParticles),
              s_(settings),
              radialGaussianBasis_(s_),
              angularBasis_(s_.angular)
    {}

    NeighborhoodExpansionCoefficientsVector<Type> expandParticle(const Particle<Type> &center,
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

        NeighborhoodExpansionCoefficientsVector<Type> singleNeighbor(1,s_);

        for (unsigned n = 1; n <= s_.radial.nmax; ++n) {
            for (unsigned l = 0; l <= s_.angular.lmax; ++l) {
                for (int m = -int(l); m <= int(l); ++m) {
                    auto coefficient = radialGaussianBasis_.computeCoefficient(n, l, centerToNeighborDistance, neighborSigma)
                                       * angularBasis_.computeCoefficient(l, m, theta, phi);

                    singleNeighbor.storeCoefficient(0, n, l, m, coefficient);
                }
            }
        }
        return singleNeighbor;
    }

    const ExpansionSettings& getSettings(){
        return s_;
    }

private:
    unsigned numberOfParticles_;
    ExpansionSettings s_;
    RadialGaussianBasis radialGaussianBasis_;
    AngularBasis angularBasis_;
};

#endif //AMOLQCPP_NEIGHBORHOODEXPANDER_H
