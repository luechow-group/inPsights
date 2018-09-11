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

        explicit Dataset(std::vector<VectorType> data)
        :
        data_(data),
        rows_(data.size()),
        cols_(data[0].size())
        {}

        std::vector<VectorType>& data() {
            return data_;
        }

        size_t rows() const {
            return rows_;
        }

        size_t cols() const {
            return cols_;
        }

    private:
        std::vector<VectorType> data_;
        size_t rows_;
        size_t cols_;
    };
}

#endif // DATASET_H
