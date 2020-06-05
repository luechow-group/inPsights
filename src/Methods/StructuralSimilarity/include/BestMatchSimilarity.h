/* Copyright (C) 2019-2020 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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

        std::deque<std::deque<PermuteeToReferenceMatch>> groupDependentMatches(
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
