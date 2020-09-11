// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "NeighborhoodExpander.h"
#include "Cutoff.h"
#include "AngularBasis.h"
#include "ParticleKit.h"
#include <SOAPSettings.h>

using namespace SOAP;

NeighborhoodExpander::NeighborhoodExpander()
        : radialGaussianBasis_(){}
        
NeighborhoodExpansion NeighborhoodExpander::expandEnvironment(const Environment& e, int expansionTypeId) const {
    NeighborhoodExpansion neighborhoodExpansion;

    auto nmax = Radial::settings.nmax();
    auto lmax = Angular::settings.lmax();
    auto sigmaAtom  = Radial::settings.sigmaAtom();
    auto radiusZero = Radial::settings.radiusZero();
    auto centerWeight = Cutoff::settings.centerWeight();

    for (const auto& neighborCoordsPair : e.selectParticles(expansionTypeId)) {

        const auto& neighborCoords = neighborCoordsPair.second;

        double weight = 1; //TODO TypeSpecific Value?
        double weightScale = Cutoff::value(neighborCoords.r);

        if (neighborCoords.r <= radiusZero)
            weight *= centerWeight;

        for (unsigned n = 1; n <= nmax; ++n) {
            for (unsigned l = 0; l <= lmax; ++l) {

                //radialGaussianBasis_.
                auto radialCoeff = radialGaussianBasis_.computeCoefficients(neighborCoords.r, sigmaAtom);
                for (int m = -int(l); m <= int(l); ++m) {

                    //TODO use TypeSpecific sigma value? Is neighbor sigma right?
                    auto coeff = radialCoeff(n-1,l)* AngularBasis::computeCoefficient(l, m, neighborCoords.theta, neighborCoords.phi)
                                 * weight*weightScale;
                    neighborhoodExpansion.storeCoefficient(n,l,m,coeff);
                }
            }
        }
    }
    return neighborhoodExpansion;
}

TypeSpecificNeighborhoodsAtOneCenter
NeighborhoodExpander::computeParticularExpansions(const Environment &e) {
    TypeSpecificNeighborhoodsAtOneCenter expansions;

    auto mode = General::settings.mode();
    if (mode == General::Mode::typeAgnostic) {
        auto noneTypeId = 0;
        expansions.emplace(noneTypeId, expandEnvironment(e, noneTypeId));
    }
    else if (mode == General::Mode::chemical || mode == General::Mode::alchemical){
        for (auto &[type, count] : ParticleKit::kit)
            expansions.emplace(type, expandEnvironment(e, type));
    } else {
        throw std::exception();
    }
    return expansions;
}

MolecularCenters
NeighborhoodExpander::computeMolecularExpansions(MolecularGeometry molecule) {
    assert(ParticleKit::isSubsetQ(molecule)
           && "The molecule must be composable from the set of particles specified in the particle kit");

    MolecularCenters expansions;

    for (unsigned k = 0; k < unsigned(molecule.numberOfEntities()); ++k) {
        auto currentEnumeratedType = molecule.findEnumeratedTypeByIndex(k);
        expansions[currentEnumeratedType] = computeParticularExpansions(Environment(molecule, currentEnumeratedType));
    }

    return expansions;
}
