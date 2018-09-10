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

        static inline double dist(const VectorType &p1, const VectorType &p2) {
            return (p1 - p2).norm();
        }

        std::shared_ptr<VantagePointTree<Scalar, dist>> get_vp() const {
            return vp_tree_;
        }

        void fit() {
            const auto &d = dset_->data();

            const double start = omp_get_wtime();

            vp_tree_ = std::make_shared<VantagePointTree<Scalar,dist >>();
            vp_tree_->create(dset_);

            const size_t dlen = d.size();

            prepareLabels(dlen);

            fit_time_ = omp_get_wtime() - start;
        }

        const std::vector<double> predictEps(size_t k) {
            const auto &d = dset_->data();

            std::vector<double> r(d.size(), 0.0);

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

        uint32_t predict(double eps, size_t min_elems){

            std::unique_ptr<std::vector<uint32_t> > candidates(new std::vector<uint32_t>());
            std::unique_ptr<std::vector<uint32_t> > new_candidates(new std::vector<uint32_t>());

            int32_t cluster_id = 0;

            std::vector<std::pair<size_t, float>> index_neigh;
            std::vector<std::pair<size_t, float>> n_neigh;

            const double start = omp_get_wtime();

            const auto &d = dset_->data();
            const size_t dlen = d.size();

            spdlog::get("console")->info("start");
            for (uint32_t pid = 0; pid < dlen; ++pid) {
                if (pid % 10000 == 0)
                    spdlog::get("console")->info("progress: pid = {0}, {1}%", pid, (float(pid) / float(dlen)) * 100);

                if (labels_[pid] >= 0)
                    continue;


                find_neighbors(d, eps, pid, index_neigh);

                spdlog::get("console")->info("Analyzing pid {0} Neigh size {1}", pid, index_neigh.size());

                if (index_neigh.size() < min_elems)
                    continue;

                labels_[pid] = cluster_id;

                spdlog::get("console")->info("pid = {0} neig = {1}", pid, index_neigh.size());

                candidates->clear();

                for (const auto &nn : index_neigh) {

                    if (labels_[nn.first] >= 0)
                        continue;

                    labels_[nn.first] = cluster_id;

                    // find_neighbors( d, eps, nn.first, n_neigh );

                    spdlog::get("console")->info("nn.first = {0} neig = {1}",nn.first, n_neigh.size());

                    // if ( n_neigh.size() >= min_elems ) {
                    candidates->push_back(nn.first);
                    // }
                }

                while (candidates->size() > 0) {
                    spdlog::get("console")->info("candidates size {}", candidates->size());
                    new_candidates->clear();

                    const float csize = float(candidates->size());

#pragma omp parallel for ordered schedule( dynamic )
                    for (size_t j = 0; j < candidates->size(); ++j) {
                        // for ( const auto& c_pid : *candidates ) {
                        std::vector<std::pair<size_t, float>> c_neigh;
                        const uint32_t c_pid = candidates->at(j);

                        spdlog::get("console")->info("c_pid = {0}, {1}", c_pid,labels_[c_pid]);

                        // if ( m_labels[c_pid] >= 0 && m_labels[c_pid] != cluster_id )
                        //     continue;

                        find_neighbors(d, eps, c_pid, c_neigh);

                        if (c_neigh.size() < min_elems)
                            continue;

                        spdlog::get("console")->info("c_pid = {0} neig = {1}", c_pid,c_neigh.size());

#pragma omp ordered
                        {
                            for (const auto &nn : c_neigh) {

                                if (labels_[nn.first] >= 0)
                                    continue;

                                labels_[nn.first] = cluster_id;

                                // find_neighbors( d, eps, nn.first, n_neigh );

                                spdlog::get("console")->info("nn.first = {0} neig = {1}", nn.first, n_neigh.size());
                                // if ( n_neigh.size() >= min_elems ) {

                                new_candidates->push_back(nn.first);
                            }
                            //if (j % 1000 == 0)
                            spdlog::get("console")->info("sub progress: j = {0} {1}% {2}", j, ( float(j)/csize )*100,new_candidates->size());
                        }
                        // }
                    }
                    spdlog::get("console")->info("new candidates = {}",new_candidates->size());

                    std::swap(candidates, new_candidates);
                }
                ++cluster_id;
            }

            predict_time_ = omp_get_wtime() - start;
            spdlog::get("console")->flush();
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
        void find_neighbors(const std::vector<VectorType> &d, double eps, uint32_t pid,
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
