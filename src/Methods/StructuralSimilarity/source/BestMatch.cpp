//
// Created by Michael Heuer on 26.09.18.
//

#include <BestMatch.h>

Eigen::PermutationMatrix<Eigen::Dynamic> BestMatch::combinePermutations(
        const Eigen::PermutationMatrix<Eigen::Dynamic> &p1,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &p2, bool flipSpinsQ) {

    long n1 = p1.size(), n2 = p2.size();
    Eigen::VectorXi combined(n1 + n2);

    if (!flipSpinsQ) {
        combined.segment(0, n1) = p1.indices().base();
        combined.segment(n1, n2) = (p2.indices().base().array() + n1);
    } else {
        combined.segment(0, n1) = (p1.indices().base().array() + n1);
        combined.segment(n1, n2) = p2.indices().base();
    }
    return Eigen::PermutationMatrix<Eigen::Dynamic>(combined);
};


BestMatch::Result BestMatch::Similarity::compare(
        const MolecularSpectrum &permutee,
        const MolecularSpectrum &reference) {
    assert(ParticleKit::isSubsetQ(permutee.molecule_) &&
           "The permutee must be a subset of the particle kit.");
    assert(ParticleKit::isSubsetQ(reference.molecule_) &&
           "The reference must be a subset of the particle kit.");

    // TODO assert that identical number of electrons and same atom geometry? Is this constraint needed? What happens with rows/cols of zero?

    assert(
            permutee.molecule_.electrons().typesVector().countOccurence(Spin::alpha) ==
            reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha)
            && "The number of alpha electrons has to match.");

    assert(
            permutee.molecule_.electrons().typesVector().countOccurence(Spin::beta) ==
            reference.molecule_.electrons().typesVector().countOccurence(Spin::beta)
            && "The number of beta electrons has to match.");

    auto nAlpha = reference.molecule_.electrons().typesVector().countOccurence(Spin::alpha);
    auto nBeta = reference.molecule_.electrons().typesVector().countOccurence(Spin::beta);

    auto N = nAlpha + nBeta;
    Eigen::MatrixXd environmentalSimilarities(N, N);

    // TODO consider identical spin flip
    TypeSpecificNeighborhoodsAtOneCenter expA, expB;
    for (unsigned i = 0; i < nAlpha; ++i) {
        EnumeratedType<int> enumeratedType_i(Spins::spinToInt(Spin::alpha), i);
        expA = permutee.molecularCenters_.find(enumeratedType_i)->second;

        for (unsigned j = 0; j < nAlpha; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::alpha), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(i, j) = LocalSimilarity::kernel(expA, expB,
                                                                      SOAPExpansion::settings.zeta());
        }
        for (unsigned j = 0; j < nBeta; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::beta), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(i, nAlpha + j) = LocalSimilarity::kernel(expA, expB,
                                                                               SOAPExpansion::settings.zeta());
        }
    }
    for (unsigned i = 0; i < nBeta; ++i) {
        EnumeratedType<int> enumeratedType_i(Spins::spinToInt(Spin::beta), i);
        expA = permutee.molecularCenters_.find(enumeratedType_i)->second;

        for (unsigned j = 0; j < nAlpha; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::alpha), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(nAlpha + i, j) = LocalSimilarity::kernel(expA, expB,
                                                                               SOAPExpansion::settings.zeta());
        }
        for (unsigned j = 0; j < nBeta; ++j) {
            EnumeratedType<int> enumeratedType_j(Spins::spinToInt(Spin::beta), j);
            expB = reference.molecularCenters_.find(enumeratedType_j)->second;
            environmentalSimilarities(nAlpha + i, nAlpha + j) = LocalSimilarity::kernel(expA, expB,
                                                                                        SOAPExpansion::settings.zeta());
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