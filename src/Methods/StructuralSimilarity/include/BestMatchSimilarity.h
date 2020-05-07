/* Copyright (C) 2019 Michael Heuer.
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


        Eigen::MatrixXd calculateEnvironmentalSimilarityMatrix(
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
                const Eigen::MatrixXd &bestMatchPermutedEnvironmentalSimilarities,
                const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch,
                double soapThreshold,
                double numericalPrecisionEpsilon = std::numeric_limits<double>::epsilon());

        double earlyExitMetric(const Eigen::MatrixXd &bestMatchPermutedEnvironmentalSimilarities);

        Eigen::MatrixXd calculateDistanceCovarianceMatrixOfSelectedIndices(
                const ElectronsVector &electronsVector,
                const std::vector<Eigen::Index> &kitSystemIndices);
    }
}

#endif //INPSIGHTS_BESTMATCHSOAPSIMILARITY_H
