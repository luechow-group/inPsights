#ifndef DENSITYBASEDSCAN_H
#define DENSITYBASEDSCAN_H

#include "VantagePointTree.h"
#include "Dataset.h"
#include <Eigen/Dense>

namespace Clustering {
    class DensityBasedScan {
    public:
        DensityBasedScan(const DensityBasedScan &) = delete;
        DensityBasedScan &operator=(const DensityBasedScan &) = delete;
        ~DensityBasedScan() = default;

        explicit DensityBasedScan(const std::shared_ptr<Dataset>& dset);

        static inline double dist(const Eigen::VectorXf &p1, const Eigen::VectorXf &p2) {
            return (p1 - p2).norm();
        }

        std::shared_ptr<VantagePointTree<Eigen::VectorXf, dist >> get_vp() const;

        void fit();

        const std::vector<double> predictEps(size_t k);

        uint32_t predict(double eps, size_t min_elems);

        void reset();

        const std::vector<int32_t> &getLabels() const;

        const double getFitTime() const;

        const double getPredictTime() const;

    private:
        void find_neighbors(const std::vector<Eigen::VectorXf> &d, double eps, uint32_t pid,
                            std::vector<std::pair<size_t, float>> &neighbors);

        void prepareLabels(size_t s);

        const std::shared_ptr<Dataset> dset_;
        std::vector<int32_t> labels_;
        std::shared_ptr<VantagePointTree<Eigen::VectorXf, dist >> vp_tree_;
        double fit_time_;
        double predict_time_;
    };

}

#endif // DENSITYBASEDSCAN_H
