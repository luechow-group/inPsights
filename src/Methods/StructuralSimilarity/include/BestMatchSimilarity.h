// Copyright (C) 2019-2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_BESTMATCHSOAPSIMILARITY_H
#define INPSIGHTS_BESTMATCHSOAPSIMILARITY_H

#include "BestMatch.h"
#include <MolecularSpectrum.h>
#include <vector>
#include <deque>

namespace BestMatch {
    namespace SOAPSimilarity {

        /*
        Eigen::MatrixXd calculateEnvironmentSimilarityMatrix(
                const SOAP::MolecularSpectrum &permutee,
                const SOAP::MolecularSpectrum &reference);*/

        DescendingMetricResult compare(
                const SOAP::MolecularSpectrum &permutee,
                const SOAP::MolecularSpectrum &reference,
                double distanceMatrixCovarianceTolerance,
                double soapThreshold,
                double numericalPrecisionEpsilon = std::numeric_limits<double>::epsilon());

        std::vector<DescendingMetricResult> getBestMatchResults(
                const SOAP::MolecularSpectrum &permutee,
                const SOAP::MolecularSpectrum &reference,
                double distanceMatrixCovarianceTolerance,
                double similarityThreshold,
                double comparisionEpsilon = std::numeric_limits<double>::epsilon());

        struct PermuteeToReferenceMatch {
            std::set<Eigen::Index> permuteeIndices;
            Eigen::Index referenceIndex;
        };

        struct GrowingPerm {
            std::set<Eigen::Index> remainingPermuteeIndices_;
            std::deque<std::pair<Eigen::Index, Eigen::Index>> chainOfSwaps_;

            GrowingPerm(std::set<Eigen::Index> remainingPermuteeIndices,
                        std::deque<std::pair<Eigen::Index, Eigen::Index>> chainOfSwaps);

            bool add(const std::pair<Eigen::Index, Eigen::Index> &envMatch);

            bool operator<(const GrowingPerm &rhs) const;
        };

        std::deque<PermuteeToReferenceMatch> findEnvironmentMatches(
                const Eigen::MatrixXd &environmentSimilarities,
                double soapThreshold,
                double numericalPrecisionEpsilon = std::numeric_limits<double>::epsilon());

        std::deque<std::deque<PermuteeToReferenceMatch>> clusterDependentMatches(
                const std::deque<PermuteeToReferenceMatch> &matches);

        std::deque<GrowingPerm> findPossiblePermutations(const std::deque<PermuteeToReferenceMatch> &dependentMatches);

        double earlyExitMetric(const Eigen::MatrixXd &bestMatchPermutedEnvironmentSimilarities);

        std::vector<Eigen::Index> permuteIndicesFromKitSystem(const std::vector<Eigen::Index> &kitSystemIndices,
                                                              const Eigen::PermutationMatrix<Eigen::Dynamic>& fromKitPermutation);

        Eigen::MatrixXd calculateDistanceCovarianceMatrixOfSelectedIndices(
                const PositionsVector &positions,
                const std::vector<Eigen::Index> &indices);
    }
}

#endif //INPSIGHTS_BESTMATCHSOAPSIMILARITY_H
