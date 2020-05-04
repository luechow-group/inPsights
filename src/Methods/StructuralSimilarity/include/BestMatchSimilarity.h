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

        Result compare(
                const SOAP::MolecularSpectrum &permutee,
                const SOAP::MolecularSpectrum &reference,
                double distanceMatrixCovarianceTolerance, double soapThreshold);

        std::vector<Result> getBestMatchResults(
                const SOAP::MolecularSpectrum &permutee,
                const SOAP::MolecularSpectrum &reference,
                double distanceMatrixCovarianceTolerance, double soapThreshold);

        std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> getBlockwiseDependentIndexPairs(
                const Eigen::MatrixXd &environmentalSimilarities,
                const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch,
                double soapThreshold);

        void varySimilarEnvironmentsInBlock(
                const MolecularGeometry &permutee,
                const MolecularGeometry &reference,
                std::deque<std::pair<Eigen::Index,Eigen::Index>> remainingIndexPairs,
                const std::deque<std::pair<Eigen::Index,Eigen::Index>>& survivingIndexPairs,
                std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> &distancePreservingEnvironmentCombinations,
                double distanceMatrixCovarianceTolerance);

        std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>> combineBlocks(
                const MolecularGeometry &permutee,
                const MolecularGeometry &reference,
                const std::deque<std::vector<std::deque<std::pair<Eigen::Index,Eigen::Index>>>> &intraBlockDistanceCombinations,
                double distanceMatrixCovarianceTolerance);

        std::vector<Eigen::Index> unblockDependentIndicesOfPreservingCombinations(
                const std::deque<std::vector<std::deque<std::pair<Eigen::Index, Eigen::Index>>>> &distancePreservingEnvironmentCombinationsOfAllBlocks);

        Eigen::MatrixXd calculateDistanceCovarianceMatrixOfSelectedIndices(
                const ElectronsVector &electronsVector,
                std::deque<Eigen::Index> kitSystemIndices);
    }
}

#endif //INPSIGHTS_BESTMATCHSOAPSIMILARITY_H
