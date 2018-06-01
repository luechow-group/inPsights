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

        if (neighborCoords.r <= ExpansionSettings::Radial::radiusZero) {
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

TypeSpecificNeighborhoodsAtOneCenter
NeighborhoodExpander::computeExpansions(const Environment &e) {

    TypeSpecificNeighborhoodsAtOneCenter expansions;

    switch (ExpansionSettings::mode) {
        case ExpansionSettings::Mode::Generic: {
            auto noneTypeId = int(GeneralStorageType::None);
            expansions.emplace(noneTypeId, expandEnvironment(e, noneTypeId));
            break;
        }
        case ExpansionSettings::Mode::TypeSpecific: {
            for(auto & type : ParticleKit::kit){
                expansions.emplace(type.first, expandEnvironment(e, type.first));
            }
            break;
        }
    }
    return expansions;
}

MolecularCenters
NeighborhoodExpander::computeExpansions(MolecularGeometry molecule) {
    assert(ParticleKit::isSubsetQ(molecule)
           && "The molecule must be composable from the set of particles specified in the particle  kit");

    MolecularCenters exp;
    
    for (int k = 0; k < molecule.atoms().numberOfEntities(); ++k) {
        auto numberedType = molecule.atoms().typesVector().getNumberedTypeByIndex(k).toIntType();
        exp[numberedType] = computeExpansions({molecule, molecule.atoms()[k].position()});
        break;
    }
    for (int k = 0; k < molecule.electrons().numberOfEntities(); ++k) {
        auto numberedType = molecule.electrons().typesVector().getNumberedTypeByIndex(k).toIntType();
        exp[numberedType] = computeExpansions({molecule, molecule.electrons()[k].position()});
        break;
    }
    return exp;
}
