//
// Created by Michael Heuer on 15.05.18.
//

#include "NeighborhoodExpander.h"
#include "CutoffFunction.h"
#include "AngularBasis.h"
#include "ParticleKit.h"
#include <iomanip>
#include <Type.h>

NeighborhoodExpander::NeighborhoodExpander()
        : radialGaussianBasis_(){}

std::complex<double>
NeighborhoodExpander::coefficient(unsigned n, unsigned l, int m, const SphericalCoordinates &coords, double weight, double weightScale,
                                  double neighborSigma) const {

    return radialGaussianBasis_.computeCoefficient(n, l, coords.r, neighborSigma)
           * AngularBasis::computeCoefficient(l, m, coords.theta, coords.phi)
           * weight*weightScale;
}

NeighborhoodExpansion NeighborhoodExpander::expandEnvironment(const Environment& e, int expansionTypeId) const {

    NeighborhoodExpansion neighborhoodExpansion;

    for (const auto& neighborCoordsPair : e.selectParticles(expansionTypeId)) {

        const auto& neighborCoords = neighborCoordsPair.second;

        double weight = 1; //TODO TypeSpecific Value? //const auto& neighbor = neighborCoordsPair.first;
        double weightScale = CutoffFunction::getWeight(neighborCoords.r);

        if (neighborCoords.r <= ZeroLimits::radiusZero) {
            weight *= ExpansionSettings::Cutoff::centerWeight;
            //TODO return something here?
        }

        for (unsigned n = 1; n <= ExpansionSettings::Radial::nmax; ++n) {
            for (unsigned l = 0; l <= ExpansionSettings::Angular::lmax; ++l) {
                for (int m = -int(l); m < int(l); ++m) {

                    //TODO use TypeSpecific sigma value?
                    auto coeff = coefficient(n, l, m, neighborCoords, weight, weightScale);//,neighborSigma);
                    neighborhoodExpansion.storeCoefficient(n,l,m,coeff);
                }
            }
        }
    }
    return neighborhoodExpansion;
}

std::map<int, NeighborhoodExpansion>
NeighborhoodExpander::computeExpansions(const Environment &e) {

    std::map<int, NeighborhoodExpansion> expansions;//TODO MODIFY

    switch (ExpansionSettings::mode) {
        case ExpansionMode::Generic: {
            auto noneTypeId = int(Type::None);
            expansions.emplace(noneTypeId, expandEnvironment(e, noneTypeId));
            break;
        }
        case ExpansionMode::TypeSpecific: {

            auto numberOfElementTypes = unsigned(ParticleKit::atomKit.size());

            for (unsigned t = 0; t < numberOfElementTypes; ++t) {
                auto typeId = int(ParticleKit::atomKit[t].first);

                expansions.emplace(typeId, expandEnvironment(e, typeId));
            };


            //TODO ALSO EXPAND W.R.T. alpha and beta electrons
            expansions.emplace(int(Spins::SpinType::alpha),
                               expandEnvironment(e, int(Spins::SpinType::alpha)));
            expansions.emplace(int(Spins::SpinType::beta),
                               expandEnvironment(e, int(Spins::SpinType::beta)));

            break;
        }
    }
    return expansions;
}