//
// Created by Michael Heuer on 15.05.18.
//

#include "NeighborhoodExpander.h"
#include "CutoffFunction.h"
#include "AngularBasis.h"
#include "ParticleKit.h"

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

        if (neighborCoords.r <= ExpansionSettings::Radial::radiusZero)
            weight *= ExpansionSettings::Cutoff::centerWeight; //TODO return something here?

        for (unsigned n = 1; n <= ExpansionSettings::Radial::nmax; ++n) {
            for (unsigned l = 0; l <= ExpansionSettings::Angular::lmax; ++l) {
                for (int m = -int(l); m <= int(l); ++m) {

                    //TODO use TypeSpecific sigma value?
                    auto coeff = coefficient(n, l, m, neighborCoords, weight, weightScale);//,neighborSigma);
                    neighborhoodExpansion.storeCoefficient(n,l,m,coeff);
                }
            }
        }
    }
    return neighborhoodExpansion;
}

TypeSpecificNeighborhoodsAtOneCenter
NeighborhoodExpander::computeParticularExpansions(const Environment &e) { // WORKS!
    TypeSpecificNeighborhoodsAtOneCenter expansions;

    switch (ExpansionSettings::mode) {
        case ExpansionSettings::Mode::generic: {
            auto noneTypeId = 0;
            expansions.emplace(noneTypeId, expandEnvironment(e, noneTypeId));
            break;
        }
        case ExpansionSettings::Mode::chemical: {
            for(auto & type : ParticleKit::kit){
                expansions.emplace(type.first, expandEnvironment(e, type.first));
            }
            break;
        }
        case ExpansionSettings::Mode::alchemical: {
            for(auto & type : ParticleKit::kit){
                expansions.emplace(type.first, expandEnvironment(e, type.first));
            }
            break;
        }
    }
    return expansions;
}

MolecularCenters
NeighborhoodExpander::computeMolecularExpansions(MolecularGeometry molecule) {
    assert(ParticleKit::isSubsetQ(molecule)
           && "The molecule must be composable from the set of particles specified in the particle kit");
    MolecularCenters exp;

    //TODO CHECK HERE FOR IDENTICAL CENTERS!
    for (unsigned k = 0; k < unsigned(molecule.numberOfEntities()); ++k) {

        // check if center was calculated already ;
        // TODO: not possible, if type specific center value is chosen which currently isn't the case;
        bool computedAlreadyQ = false;

        NumberedType<int> existingNumberedType;
        for (unsigned i = 0; i < k; ++i) {
            if((molecule[i].position()-molecule[k].position()).norm() <= ExpansionSettings::Radial::radiusZero){
                existingNumberedType = molecule.findNumberedTypeByIndex(i);
                computedAlreadyQ = true;
                //std::cout << "found " << existingNumberedType << std::endl;
                break;
            }
        }
        auto newNumberedType = molecule.findNumberedTypeByIndex(k);
        if(computedAlreadyQ){
            exp[newNumberedType] = exp[existingNumberedType];
        } else {
            exp[newNumberedType] = computeParticularExpansions(Environment(molecule, molecule[k].position()));
        }
    }

    return exp;
}
