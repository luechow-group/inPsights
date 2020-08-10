/*
 * Copyright (c) 2017, Vladimir Skvortsov
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the FreeBSD Project.
 */
/*
 * Copyright (C) 2018 Michael Heuer.
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
