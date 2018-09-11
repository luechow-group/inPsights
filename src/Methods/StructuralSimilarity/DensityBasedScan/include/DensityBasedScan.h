#ifndef DENSITYBASEDSCAN_H
#define DENSITYBASEDSCAN_H

#include "VantagePointTree.h"

namespace Clustering {
    template <typename Scalar>
    class DensityBasedScan {
        using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    public:
        DensityBasedScan(const DensityBasedScan &) = delete;
        DensityBasedScan &operator=(const DensityBasedScan &) = delete;
        ~DensityBasedScan() = default;

        explicit DensityBasedScan(const std::vector<VectorType>& data, Scalar similarityDistance = 1e-7)
        :
        data_(data),
        vpTree_(data, similarityDistance),
        labels_(data.size(),-1) {};

        static inline Scalar dist(const VectorType &p1, const VectorType &p2) {
            return (p1 - p2).norm();
        }

        const std::vector<Scalar> predictEps(size_t k) {
            std::vector<Scalar> r(data_.size(), 0.0);

            omp_set_dynamic(1);

#pragma omp parallel for
            for (size_t i = 0; i < data_.size(); ++i) {
                std::vector<std::pair<size_t, Scalar>> nlist;

                vpTree_.searchByK(data_[i], k, nlist, true);

                if (nlist.size() >= k) {
                    r[i] = nlist[0].second;
                }
            }

            std::sort(r.begin(), r.end());

            return std::move(r);
        }

        int32_t predict(Scalar eps, size_t minPts){
            std::vector<Eigen::Index> candidates, newCandidates;
            std::vector<std::pair<size_t, Scalar>> index_neigh, n_neigh;

            const auto start = omp_get_wtime();
            const size_t dlen = data_.size();

            int32_t clusterId = 0;
            for (size_t pid = 0; pid < dlen; ++pid) {

                if (labels_[pid] >= 0) continue;

                findNeighbors(eps, pid, index_neigh);

                if (index_neigh.size() < minPts) continue;

                labels_[pid] = clusterId;

                candidates.clear();

                for (const auto &nn : index_neigh) {

                    if (labels_[nn.first] >= 0) continue;

                    labels_[nn.first] = clusterId;

                    candidates.push_back(nn.first);
                }

                while (candidates.size() > 0) {
                    newCandidates.clear();

                    const Scalar csize = Scalar(candidates.size());

#pragma omp parallel for ordered schedule( dynamic )
                    for (size_t j = 0; j < candidates.size(); ++j) {
                        std::vector<std::pair<size_t, Scalar>> c_neigh;
                        const Eigen::Index c_pid = candidates.at(j);

                        findNeighbors(eps, c_pid, c_neigh);

                        if (c_neigh.size() < minPts) continue;

#pragma omp ordered
                        {
                            for (const auto &nn : c_neigh) {

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

        void reset() {
            vpTree_.reset();
            labels_.clear();
        }

        const std::vector<int32_t> &getLabels() const {
            return labels_;
        }

    private:
        void findNeighbors(Scalar eps, Eigen::Index pid, std::vector<std::pair<size_t, Scalar>> &neighbors) {
            neighbors.clear();
            vpTree_.searchByDistance(data_[pid], eps, neighbors);
        }

        const std::vector<VectorType>& data_;
        VantagePointTree<Scalar, dist > vpTree_;
        std::vector<int32_t> labels_;
    };

}

#endif // DENSITYBASEDSCAN_H
