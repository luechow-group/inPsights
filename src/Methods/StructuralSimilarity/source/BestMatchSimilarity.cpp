//
// Created by heuer on 03.04.19.
//

#include "BestMatchSimilarity.h"
#include <LocalSimilarity.h>
#include <Hungarian.h>
#include <BestMatchSimilarity.h>

using namespace SOAP;

BestMatch::Result BestMatch::Similarity::compare(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference) {

    assert(ParticleKit::isSubsetQ(permutee.molecule_)
           && "The permutee must be a subset of the particle kit.");
    assert(ParticleKit::isSubsetQ(reference.molecule_)
           && "The reference must be a subset of the particle kit.");

    // TODO assert that identical number of electrons and same atom geometry? Is this constraint needed? What happens with rows/cols of zero?

    assert(permutee.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
           == reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
           && "The number of alpha electrons has to match.");

    assert(permutee.molecule_.electrons().typesVector().countOccurence(Spin::beta)
           == reference.molecule_.electrons().typesVector().countOccurence(Spin::beta)
           && "The number of beta electrons has to match.");

    auto nAlpha = reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha);
    auto nBeta = reference.molecule_.electrons().typesVector().countOccurence(Spin::beta);

    auto N = nAlpha + nBeta;
    Eigen::MatrixXd environmentalSimilarities(N, N);

    // TODO consider identical spin flip?

    TypeSpecificNeighborhoodsAtOneCenter expA, expB;
    for (unsigned i = 0; i < nAlpha; ++i) {
        EnumeratedType<int> enumeratedType_i(Spins::spinToInt(Spin::alpha), i);
        expA = permutee.molecularCenters_.find(enumeratedType_i)->second;

        for (unsigned j = 0; j < nAlpha; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::alpha), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(i, j) = LocalSimilarity::kernel(expA, expB, General::settings.zeta());
        }
        for (unsigned j = 0; j < nBeta; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::beta), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(i, nAlpha + j) = LocalSimilarity::kernel(expA, expB, General::settings.zeta());
        }
    }
    for (unsigned i = 0; i < nBeta; ++i) {
        EnumeratedType<int> enumeratedType_i(Spins::spinToInt(Spin::beta), i);
        expA = permutee.molecularCenters_.find(enumeratedType_i)->second;

        for (unsigned j = 0; j < nAlpha; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::alpha), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(nAlpha + i, j) = LocalSimilarity::kernel(expA, expB, General::settings.zeta());
        }
        for (unsigned j = 0; j < nBeta; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::beta), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(nAlpha + i, nAlpha + j) = LocalSimilarity::kernel(expA, expB,
                    General::settings.zeta());
        }
    }

    Eigen::PermutationMatrix<Eigen::Dynamic> bestMatch = Hungarian<double>::findMatching(
            environmentalSimilarities, Matchtype::MAX);

    // best-match permute columns and sum diagonal elements
    double simMetric = (environmentalSimilarities * bestMatch).diagonal().sum() / N;

    //restore the original order before the particle kit permutations
    auto permuteeToKit = ParticleKit::toKitPermutation(permutee.molecule_.electrons());
    auto referenceFromKit = ParticleKit::fromKitPermutation(reference.molecule_.electrons());

    return {simMetric, referenceFromKit * bestMatch * permuteeToKit};
}

BestMatch::Result BestMatch::Similarity::compare(
        MolecularGeometry permutee, const MolecularGeometry &reference,
        bool spinSpecificQ, bool flipSpinsQ) {

    ParticleKit::create(reference);
    MolecularSpectrum permuteeSpectrum, referenceSpectrum;


    if (spinSpecificQ) {
        ParticleKit::create(reference);
        referenceSpectrum = MolecularSpectrum(reference);

        if (flipSpinsQ)
            permutee.electrons().typesVector().flipSpins();

        assert(ParticleKit::isSubsetQ(permutee) && "The permutee must be a subset of the particle kit.");
        permuteeSpectrum = MolecularSpectrum(permutee);

    } else {
        // slower variant
        /*SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
        // make alpha and beta spins identical
        SOAPExpansion::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;*/

        // Trick: Instead of using SOAPExpansion::Mode::alchemical construct a new MolecularGeometry with an
        // ElectronsVector where all spins are alpha electrons
        auto permuteeUnspecific = MolecularGeometry(permutee.atoms(), {permutee.electrons().positionsVector()});
        auto referenceUnspecific = MolecularGeometry(reference.atoms(), {reference.electrons().positionsVector()});

        ParticleKit::create(referenceUnspecific);
        referenceSpectrum = MolecularSpectrum(referenceUnspecific);

        assert(ParticleKit::isSubsetQ(permuteeUnspecific) && "The permutee must be a subset of the particle kit.");
        permuteeSpectrum = MolecularSpectrum(permuteeUnspecific);
    }

    return compare(permuteeSpectrum, referenceSpectrum);
}
