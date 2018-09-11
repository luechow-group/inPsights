#ifndef DATASET_H
#define DATASET_H

#include <unordered_map>
#include <vector>
#include <memory>

#include <Eigen/Core>
#include <fstream>
#include <spdlog/spdlog.h>

namespace Clustering {
    template <typename Scalar>
    class Dataset {
        using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;
    public:

        Dataset(const Dataset&) = delete;
        Dataset& operator=(const Dataset&) = delete;

        static std::shared_ptr<Dataset> create() {
            return std::make_shared<Dataset>();
        }

        Dataset()
        : rows_(0), cols_(0) {}

        Dataset(std::vector<VectorType> data)
        :
        data_(data),
        labels_(data.size(),0),
        rows_(data.size()),
        cols_(data[0].size())
        {}

        std::vector<VectorType>& data() {
            return data_;
        }

        const std::string get_label(size_t id) const{
            auto r = reverse_labels_.find(labels_[id]);
            if (r == reverse_labels_.end()) {
                return "Unknown";
            }
            return r->second;
        }

        size_t num_points() const {
            return cols_ * rows_;
        }

        size_t rows() const {
            return rows_;
        }

        size_t cols() const {
            return cols_;
        }

    protected:
        std::vector<VectorType> data_;
        std::vector<size_t> labels_;
        size_t rows_;
        size_t cols_;
        std::unordered_map<size_t, std::string> reverse_labels_;
    };
}

#endif // DATASET_H
