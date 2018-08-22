#ifndef DATASET_H
#define DATASET_H

#include <unordered_map>
#include <vector>
#include <memory>

#include <Eigen/Core>


namespace Clustering {
    class Dataset {
    public:

        Dataset(const Dataset&) = delete;
        Dataset& operator=(const Dataset&) = delete;

        static std::shared_ptr<Dataset> create();

        Dataset();

        //void gen_cluster_data(size_t features_num, size_t elements_num);

        bool load_csv(const std::string &csv_file_path);

        std::vector<Eigen::VectorXf>& data();

        const std::string get_label(size_t id) const;

        size_t num_points() const;

        size_t rows() const;

        size_t cols() const;

    protected:
        std::vector<Eigen::VectorXf> data_;
        std::vector<size_t> labels_;
        size_t rows_;
        size_t cols_;
        std::unordered_map<size_t, std::string> reverse_labels_;
    };
}

#endif // DATASET_H
