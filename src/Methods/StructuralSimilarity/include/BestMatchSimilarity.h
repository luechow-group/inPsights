//
// Created by heuer on 03.04.19.
//

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
                double similarityRadius, double soapThreshold);

        std::vector<Result> getBestMatchResults(
                const SOAP::MolecularSpectrum &permutee,
                const SOAP::MolecularSpectrum &reference,
                double similarityRadius, double soapThreshold);

        std::vector<std::deque<Eigen::Index>> getListOfDependentIndicesLists(
                const Eigen::MatrixXd &environmentalSimilarities, double soapThreshold);

        void varySimilarEnvironments(
                const MolecularGeometry &permutee,
                const MolecularGeometry &reference,
                std::deque<Eigen::Index> dependentIndices,
                std::deque<Eigen::Index> surviving,
                std::vector<std::deque<Eigen::Index>> &allPerms,
                double similarityRadius);

        std::vector<std::deque<Eigen::Index>> combineBlocks(
                const MolecularGeometry &permutee,
                const MolecularGeometry &reference,
                const std::deque<std::vector<std::deque<Eigen::Index>>> &distancePreservingEnvironmentCombinationsOfRemainingBlocks,
                double similarityRadius);

        std::vector<Eigen::Index> obtainIndexReorderingPermutationOverAllBlocks(
                const std::deque<std::vector<std::deque<Eigen::Index>>> &distancePreservingEnvironmentCombinationsOfRemainingBlocks);

        Eigen::MatrixXd indicesBlockCovariance(
                const ElectronsVector &electronsVector,
                std::deque<Eigen::Index> indices);
    }
}

#endif //INPSIGHTS_BESTMATCHSOAPSIMILARITY_H
