// Copyright (C) 2017 Vladimir Skvortsov
// SPDX-License-Identifier: BSD-2-Clause-Views

// Copyright (C) 2018 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef DENSITYBASEDSCAN_H
#define DENSITYBASEDSCAN_H

#include "VantagePointTree.h"
#include <BestMatch.h>



struct ClusterLabels{
    int32_t numberOfClusters;
    std::vector<int32_t> labels;
};


template<typename Scalar,
        typename VectorType,
        Scalar (*distance)(const VectorType &,const VectorType &)>
class DensityBasedScan {
public:
    DensityBasedScan(const DensityBasedScan &) = delete;

    DensityBasedScan &operator=(const DensityBasedScan &) = delete;

    ~DensityBasedScan() = default;

    explicit DensityBasedScan(const std::vector<VectorType> &data, Scalar similarityDistance = 1e-7)
            :
            data_(data),
            vpTree_(data, similarityDistance) {};

    // predict eps for a fixed number of clusters
    std::vector<Scalar> predictEps(size_t k) {
        std::vector<Scalar> predictedEps(data_.size(), 0.0);

        omp_set_dynamic(1);

#pragma omp parallel shared(k, predictedEps, vpTree_)
#pragma omp for
        for (size_t i = 0; i < data_.size(); ++i) {
            std::vector<std::pair<size_t, Scalar>> neighborList;

            vpTree_.searchByK(data_[i], k, neighborList, true);

            if (neighborList.size() >= k) {
                predictedEps[i] = neighborList[0].second;
            }
        }

        std::sort(predictedEps.begin(), predictedEps.end());

        return predictedEps;
    }

    ClusterLabels findClusters(Scalar eps, size_t minPts) {
        std::vector<Eigen::Index> candidates, newCandidates;
        std::vector<std::pair<size_t, Scalar>> neighborIndices;
        std::vector<int32_t> labels(data_.size(), -1);
        
        int32_t clusterId = 0;
        for (size_t pointIndex = 0; pointIndex < data_.size(); ++pointIndex) {

            if (labels[pointIndex] >= 0) continue;

            findNeighbors(eps, pointIndex, neighborIndices);

            if (neighborIndices.size() < minPts) continue;

            labels[pointIndex] = clusterId;

            candidates.clear();

            for (const auto &nn : neighborIndices) {

                if (labels[nn.first] >= 0) continue;

                labels[nn.first] = clusterId;
                candidates.push_back(nn.first);
            }

            while (!candidates.empty()) {
                newCandidates.clear();

#pragma omp parallel for ordered schedule( dynamic )
                for (size_t j = 0; j < candidates.size(); ++j) {
                    std::vector<std::pair<size_t, Scalar>> candidateNeighbors;
                    const Eigen::Index candidateIndex = candidates.at(j);

                    findNeighbors(eps, candidateIndex, candidateNeighbors);

                    if (candidateNeighbors.size() < minPts) continue;
#pragma omp ordered
                    {
                        for (const auto &nn : candidateNeighbors) {
                            if (labels[nn.first] >= 0) continue;
                            labels[nn.first] = clusterId;
                            newCandidates.push_back(nn.first);
                        }
                    }
                }

                std::swap(candidates, newCandidates);
            }
            ++clusterId;
        }

        auto nClusters = clusterId;
        
        return {nClusters, labels};
    }

private:
    void findNeighbors(Scalar eps, Eigen::Index i, std::vector<std::pair<size_t, Scalar>> &neighbors) {
        neighbors.clear();
        vpTree_.searchByDistance(data_[i], eps, neighbors);
    }

    const std::vector<VectorType> &data_;
    VantagePointTree<Scalar,VectorType,distance> vpTree_;
};

#endif // DENSITYBASEDSCAN_H
