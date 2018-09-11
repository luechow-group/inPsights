#ifndef DENSITYBASEDSCAN_H
#define DENSITYBASEDSCAN_H

#include "VantagePointTree.h"
#include "Dataset.h"
#include <Eigen/Dense>

namespace Clustering {
    template <typename Scalar>
    class DensityBasedScan {
        using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    public:
        DensityBasedScan(const DensityBasedScan &) = delete;
        DensityBasedScan &operator=(const DensityBasedScan &) = delete;
        ~DensityBasedScan() = default;

        explicit DensityBasedScan(const std::shared_ptr<Dataset<Scalar>>& dset)
        :
        dset_(dset),
        fit_time_(.0),
        predict_time_(.0)
        {};

        static inline Scalar dist(const VectorType &p1, const VectorType &p2) {
            return (p1 - p2).norm();
        }

        std::shared_ptr<VantagePointTree<Scalar, dist>> get_vp() const {
            return vp_tree_;
        }

        void fit() {
            const auto &d = dset_->data();

            const auto start = omp_get_wtime();

            vp_tree_ = std::make_shared<VantagePointTree<Scalar,dist >>();
            vp_tree_->create(dset_);

            const size_t dlen = d.size();

            prepareLabels(dlen);

            fit_time_ = omp_get_wtime() - start;
        }

        const std::vector<Scalar> predictEps(size_t k) {
            const auto &d = dset_->data();

            std::vector<Scalar> r(d.size(), 0.0);

            omp_set_dynamic(1);

#pragma omp parallel for
            for (size_t i = 0; i < d.size(); ++i) {
                std::vector<std::pair<size_t, Scalar>> nlist;

                vp_tree_->searchByK(d[i], k, nlist, true);

                if (nlist.size() >= k) {
                    r[i] = nlist[0].second;
                }
            }

            std::sort(r.begin(), r.end());

            return std::move(r);
        }

        Eigen::Index predict(Scalar eps, size_t min_elems){

            std::unique_ptr<std::vector<Eigen::Index> > candidates(new std::vector<Eigen::Index>());
            std::unique_ptr<std::vector<Eigen::Index> > new_candidates(new std::vector<Eigen::Index>());

            int32_t cluster_id = 0;

            std::vector<std::pair<size_t, Scalar>> index_neigh;
            std::vector<std::pair<size_t, Scalar>> n_neigh;

            const auto start = omp_get_wtime();
            const auto &d = dset_->data();
            const size_t dlen = d.size();

            for (size_t pid = 0; pid < dlen; ++pid) {

                if (labels_[pid] >= 0)
                    continue;

                find_neighbors(d, eps, pid, index_neigh);

                if (index_neigh.size() < min_elems) continue;

                labels_[pid] = cluster_id;

                candidates->clear();

                for (const auto &nn : index_neigh) {

                    if (labels_[nn.first] >= 0) continue;

                    labels_[nn.first] = cluster_id;

                    candidates->push_back(nn.first);
                }

                while (candidates->size() > 0) {
                    new_candidates->clear();

                    const Scalar csize = Scalar(candidates->size());

#pragma omp parallel for ordered schedule( dynamic )
                    for (size_t j = 0; j < candidates->size(); ++j) {
                        std::vector<std::pair<size_t, Scalar>> c_neigh;
                        const Eigen::Index c_pid = candidates->at(j);

                        find_neighbors(d, eps, c_pid, c_neigh);

                        if (c_neigh.size() < min_elems) continue;

#pragma omp ordered
                        {
                            for (const auto &nn : c_neigh) {

                                if (labels_[nn.first] >= 0) continue;

                                labels_[nn.first] = cluster_id;

                                new_candidates->push_back(nn.first);
                            }
                        }
                    }

                    std::swap(candidates, new_candidates);
                }
                ++cluster_id;
            }

            predict_time_ = omp_get_wtime() - start;

            return cluster_id;
        }

        void reset() {
            vp_tree_.reset();
            labels_.clear();
        }

        const std::vector<int32_t> &getLabels() const {
            return labels_;
        }

        const double getFitTime() const {
            return fit_time_;
        }

        const double getPredictTime() const {
            return predict_time_;
        }

    private:
        void find_neighbors(const std::vector<VectorType> &d, Scalar eps, Eigen::Index pid,
                            std::vector<std::pair<size_t, Scalar>> &neighbors){
            neighbors.clear();
            vp_tree_->searchByDistance(d[pid], eps, neighbors);
        }

        void prepareLabels(size_t s) {
            labels_.resize(s);

            for (auto &l : labels_) {
                l = -1;
            }
        }

        const std::shared_ptr<Dataset<Scalar>> dset_;
        std::vector<int32_t> labels_;
        std::shared_ptr<VantagePointTree<Scalar, dist >> vp_tree_;
        double fit_time_;
        double predict_time_;
    };

}

#endif // DENSITYBASEDSCAN_H
