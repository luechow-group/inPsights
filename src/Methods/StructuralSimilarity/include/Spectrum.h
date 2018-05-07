//
// Created by Michael Heuer on 18.04.18.
//

#ifndef AMOLQCPP_SPECTRUM_H
#define AMOLQCPP_SPECTRUM_H

#include "Cutoff.h"
#include "SoapExpansion.h"
#include <ParticlesVector.h>
#include <vector>

template <typename Type>
class Spectrum{
public:

    explicit Spectrum(const ParticlesVector<Type>& particlesVector,
                      const ExpansionSettings& settings = ExpansionSettings::defaults())
            : particlesVector_(particlesVector),
              coefficientsVector_(unsigned(particlesVector_.numberOfEntities()), settings),
              soapExpansion_(settings)
    {}

    void compute(){
        for (unsigned i = 0; i < particlesVector_.numberOfEntities(); ++i) {
            expandParticlesInNeighborhood(i);
        }
    }

    void expandParticlesInNeighborhood(unsigned centerParticleId){
        assert(centerParticleId >= 0
               && "The center particle ID must be greater or equal to zero.");
        assert(centerParticleId < particlesVector_.numberOfEntities()
               && "The center particle ID must be smaller than the total number of particles.");

        for (unsigned j = 0; j < particlesVector_.numberOfEntities(); ++j) {

            /*//TODO arbitrary? change this */ double neighborSigma = soapExpansion_.getSettings().radial.sigmaAtom;

            double centerToNeighborDistance = Cutoff::distance(particlesVector_[j].position(),
                                                               particlesVector_[centerParticleId].position());

            // skip this iteration if particle i is outside the cutoff radius
            if (!cutoffFunction_.withinCutoffRadiusQ(centerToNeighborDistance)){
                continue;
            } else {
                double weight = 1;
                double weightScale = cutoffFunction_.getWeight(centerToNeighborDistance);

                if (j == centerParticleId) weight *= cutoffFunction_.getCenterWeight();

                auto atomicCoefficientsVector = soapExpansion_.expandParticle(particlesVector_[centerParticleId],
                                                                              particlesVector_[j], neighborSigma);
                atomicCoefficientsVector *= weight*weightScale;

                coefficientsVector_.storeParticleExpansion(j, atomicCoefficientsVector);
            }
        }
    }

    CoefficientsVector<Type> getCoefficients(){
        return coefficientsVector_;
    }


private:
    ParticlesVector<Type> particlesVector_;
    SoapExpansion<Type> soapExpansion_;
    CoefficientsVector<Type> coefficientsVector_;
    Cutoff cutoffFunction_;
};

#endif //AMOLQCPP_SPECTRUM_H
