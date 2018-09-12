#ifndef DENSITYBASEDSCAN_H
#define DENSITYBASEDSCAN_H

#include "VantagePointTree.h"

template<typename Scalar>
class DensityBasedScan {
    using VectorType = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
public:
    DensityBasedScan(const DensityBasedScan &) = delete;

    DensityBasedScan &operator=(const DensityBasedScan &) = delete;

    ~DensityBasedScan() = default;

    explicit DensityBasedScan(const std::vector<VectorType> &data, Scalar similarityDistance = 1e-7)
            :
            data_(data),
            vpTree_(data, similarityDistance),
            labels_(data.size(), -1) {};

    static inline Scalar euclideanDistance(const VectorType &p1, const VectorType &p2) {
        return (p1 - p2).norm();
    }

    // predict eps for a fixed number of clusters
    const std::vector<Scalar> predictEps(size_t k) {
        std::vector<Scalar> predictedEps(data_.size(), 0.0);

        omp_set_dynamic(1);

#pragma omp parallel for
        for (size_t i = 0; i < data_.size(); ++i) {
            std::vector<std::pair<size_t, Scalar>> neighborList;

            vpTree_.searchByK(data_[i], k, neighborList, true);

            if (neighborList.size() >= k) {
                predictedEps[i] = neighborList[0].second;
            }
        }

        std::sort(predictedEps.begin(), predictedEps.end());

        return std::move(predictedEps);
    }

    int32_t findClusters(Scalar eps, size_t minPts) {
        std::vector<Eigen::Index> candidates, newCandidates;
        std::vector<std::pair<size_t, Scalar>> neighborIndices, n_neigh;

        const auto start = omp_get_wtime();

        int32_t clusterId = 0;
        for (size_t pointIndex = 0; pointIndex < data_.size(); ++pointIndex) {

            if (labels_[pointIndex] >= 0) continue;

            findNeighbors(eps, pointIndex, neighborIndices);

            if (neighborIndices.size() < minPts) continue;

            labels_[pointIndex] = clusterId;

            candidates.clear();

            for (const auto &nn : neighborIndices) {

                if (labels_[nn.first] >= 0) continue;

                labels_[nn.first] = clusterId;
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
                            if (labels_[nn.first] >= 0) continue;
                            labels_[nn.first] = clusterId;
                            newCandidates.push_back(nn.first);
                        }
                    }
                }

                std::swap(candidates, newCandidates);
            }
            ++clusterId;
        }

        return clusterId;
    }

    const std::vector<int32_t> &getLabels() const {
        return labels_;
    }

private:
    void findNeighbors(Scalar eps, Eigen::Index i, std::vector<std::pair<size_t, Scalar>> &neighbors) {
        neighbors.clear();
        vpTree_.searchByDistance(data_[i], eps, neighbors);
    }

    const std::vector<VectorType> &data_;
    VantagePointTree<Scalar, euclideanDistance> vpTree_;
    std::vector<int32_t> labels_;
};

#endif // DENSITYBASEDSCAN_H
