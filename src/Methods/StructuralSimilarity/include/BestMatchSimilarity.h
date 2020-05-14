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

        Eigen::MatrixXd calculateEnvironmentSimilarityMatrix(
                const SOAP::MolecularSpectrum &permutee,
                const SOAP::MolecularSpectrum &reference);

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

        std::vector<std::deque<std::pair<Eigen::Index, Eigen::Index>>> findEquivalentEnvironments(
                const Eigen::MatrixXd &bestMatchPermutedEnvironmentSimilarities,
                const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch,
                double soapThreshold,
                double numericalPrecisionEpsilon = std::numeric_limits<double>::epsilon());

        struct PermuteeEnvsToReferenceEnvMatch {
            std::set<Eigen::Index> permuteeEnvsIndices;
            Eigen::Index referenceEnvIndex;
        };

        struct GrowingPerm {
            std::set<Eigen::Index> remainingPermuteeIndices_;
            std::deque<std::pair<Eigen::Index, Eigen::Index>> chainOfSwaps_;

            GrowingPerm(const std::set<Eigen::Index> &remainingPermuteeIndices,
                        const std::deque<std::pair<Eigen::Index, Eigen::Index>> &chainOfSwaps);

            bool add(const std::pair<Eigen::Index, Eigen::Index> &envMatch);
        };

        std::deque<PermuteeEnvsToReferenceEnvMatch> findEnvironmentMatches(
                const Eigen::MatrixXd &environmentSimilarities,
                double soapThreshold,
                double numericalPrecisionEpsilon = std::numeric_limits<double>::epsilon());

        std::deque<std::deque<PermuteeEnvsToReferenceEnvMatch>> groupDependentMatches(
                const std::deque<PermuteeEnvsToReferenceEnvMatch> &matches);

        std::deque<GrowingPerm>
        findPossiblePermutations(const std::deque<PermuteeEnvsToReferenceEnvMatch> &dependentMatches);

        double earlyExitMetric(const Eigen::MatrixXd &bestMatchPermutedEnvironmentSimilarities);

        Eigen::MatrixXd calculateDistanceCovarianceMatrixOfSelectedIndices(
                const ElectronsVector &electronsVector,
                const std::vector<Eigen::Index> &kitSystemIndices);
    }
}

#endif //INPSIGHTS_BESTMATCHSOAPSIMILARITY_H
