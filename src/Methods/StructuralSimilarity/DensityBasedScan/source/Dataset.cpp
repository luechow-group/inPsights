//
// Created by Michael Heuer on 18.08.18.
//

#include "Dataset.h"

#include <fstream>
#include <spdlog/spdlog.h>

namespace Clustering {

    std::shared_ptr<Dataset> Dataset::create() {
        return std::make_shared<Dataset>();
    }

    Dataset::Dataset()
            : rows_(0), cols_(0) {
    }

    /*void Dataset::gen_cluster_data(size_t features_num, size_t elements_num) {
        _data.clear();
        _labels.resize(elements_num);

        for (size_t i = 0; i < elements_num; ++i) {
            Eigen::VectorXf col_vector(features_num);
            for (size_t j = 0; j < features_num; ++j) {
                col_vector(j) = (-1.0 + rand() * (2.0) / RAND_MAX);
            }
            _data.emplace_back(col_vector);
        }
    }*/

    bool Dataset::load_csv(const std::string &csv_file_path) {
        std::ifstream in(csv_file_path);

        if (!in.is_open()) {
            spdlog::get("console")->error("File not opened");
            return false;
        }

        std::unordered_map<std::string, size_t> known_labels;
        reverse_labels_.clear();

        std::string line;

        std::vector<float> rowCache;
        labels_.clear();

        size_t label_idx = 0;

        while (std::getline(in, line)) {
            if (line.empty()) { //Orig if( !line.size() ) {
                continue;
            }

            rowCache.clear();

            const char *ptr = line.c_str();
            size_t len = line.length();

            const char *start = ptr;
            for (size_t i = 0; i < len; ++i) {

                if (ptr[i] == ',') {
                    rowCache.push_back(std::atof(start));
                    start = ptr + i + 1;
                }
            }

            const std::string label_str(start);

            auto r = known_labels.find(start);

            size_t found_label = label_idx;

            if (r == known_labels.end()) {
                known_labels.insert(std::make_pair(label_str, label_idx));
                reverse_labels_.insert(std::make_pair(label_idx, label_str));

                spdlog::get("console")->info("Found new label {0}",label_str);

                ++label_idx;
            } else {
                found_label = r->second;
            }

            if (!cols_) {
                cols_ = rowCache.size();
            } else {

                if (cols_ != rowCache.size()) {
                    spdlog::get("console")->error("Corrupted line \"{0}\"", line);
                    spdlog::get("console")->error("Row size = {0} line size = {1}", cols_, rowCache.size());

                    continue;
                }
            }

            Eigen::VectorXf col_vector(rowCache.size());
            for (size_t i = 0; i < rowCache.size(); ++i) {
                col_vector(i) = rowCache[i];
            }

            data_.emplace_back(col_vector);
            labels_.push_back(found_label);

            ++rows_;
        }

        in.close();

        assert(data_.size() == labels_.size());

        return true;
    }

    std::vector<Eigen::VectorXf>& Dataset::data() {
        return data_;
    }

    const std::string Dataset::get_label(size_t id) const {
        auto r = reverse_labels_.find(labels_[id]);
        if (r == reverse_labels_.end()) {
            return "Unknown";
        }
        return r->second;
    }

    size_t Dataset::num_points() const {
        return cols_ * rows_;
    }

    size_t Dataset::rows() const {
        return rows_;
    }

    size_t Dataset::cols() const {
        return cols_;
    }

};
