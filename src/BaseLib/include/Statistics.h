//
// Created by Michael Heuer on 14.09.18.
//

#ifndef INPSIGHTS_STATISTICS_H
#define INPSIGHTS_STATISTICS_H

#include <Eigen/Core>
#include <yaml-cpp/yaml.h>
#include <iostream>

namespace Statistics {

    template<typename Derived>
    class MeanErrorPair {
    public:
        MeanErrorPair(Derived mean, Derived error)
                : mean_(std::move(mean)), error_(std::move(error)) {};

        Derived mean_, error_;
    };

    template<typename Derived, typename WeightType = unsigned>
    class RunningStatistics {

    public:
        RunningStatistics()
                :
                initializedQ_(false), squaredTotalWeightInitializedQ_(false),
                totalWeight_(0), squaredTotalWeight_(0),
                mean_(), lastMean_(),
                unnormalisedVariance_(), lastUnnormalisedVariance_(),
                cwiseMin_(), cwiseMax_() {}

        RunningStatistics(
                Derived mean, Derived standardError,
                Derived cwiseMin, Derived cwiseMax,
                WeightType totalWeight, WeightType squaredTotalWeight = 0)
                :
                initializedQ_(true),
                squaredTotalWeightInitializedQ_(squaredTotalWeight == 0),
                totalWeight_(totalWeight),
                squaredTotalWeight_(squaredTotalWeightInitializedQ_? squaredTotalWeight : totalWeight*totalWeight),
                mean_(),
                lastMean_(std::move(mean)),
                unnormalisedVariance_(),
                lastUnnormalisedVariance_(calculateVarianceFromStandardError(standardError, totalWeight)),
                cwiseMin_(cwiseMin),
                cwiseMax_(cwiseMax) {}

        RunningStatistics(Derived mean, Derived standardError, WeightType totalWeight, WeightType squaredTotalWeight = 0)
                : RunningStatistics(mean, standardError,
                                    mean - 3 * standardError, mean + 3 * standardError,
                                    totalWeight, squaredTotalWeight) {}


        void reset() {
            initializedQ_ = false;
            squaredTotalWeightInitializedQ_ = false;
            totalWeight_ = 0;
            squaredTotalWeight_ = 0;
            mean_.setZero();
            lastMean_.setZero();
            unnormalisedVariance_.setZero();
            lastUnnormalisedVariance_.setZero();
            cwiseMin_.setZero();
            cwiseMax_.setZero();
        }

        void add(const Derived &sample, WeightType w = 1) {
            totalWeight_ += w;
            squaredTotalWeight_ += w * w;

            if (!initializedQ_) {
                initialize(sample);
            } else {
                Derived delta = sample - lastMean_;
                mean_ = lastMean_ + (delta * w / totalWeight_).matrix();
                unnormalisedVariance_ =
                        lastUnnormalisedVariance_ + w * (delta.array() * (sample - mean_).array()).matrix();

                lastMean_ = mean_;
                lastUnnormalisedVariance_ = unnormalisedVariance_;
                cwiseMin_ = cwiseMin_.cwiseMin(sample);
                cwiseMax_ = cwiseMax_.cwiseMax(sample);
            }
        }

        Derived mean() const {
            return mean_;
        }

        Derived cwiseMin() const {
            return cwiseMin_;
        }

        Derived cwiseMax() const {
            return cwiseMax_;
        }

        Derived variance() const {
            return unbiasedSampleVariance();
        }

        Derived standardDeviation() const {
            return unbiasedSampleStandardDeviation();
        }

        Derived standardError() const {
            return standardDeviation() / std::sqrt(getTotalWeight());
        }

        WeightType getTotalWeight() const {
            return totalWeight_;
        }

        void toYaml(YAML::Emitter &out, bool excludeSelfinteractionQ = false) const {
            using namespace YAML;

            out << BeginMap;
            for (Eigen::Index i = 0; i < mean_.rows() - (excludeSelfinteractionQ ? 1 : 0); ++i) {
                out << Key << i << Value;

                if (mean_.cols() > 1) { // Matrix
                    out << BeginMap;
                    for (Eigen::Index j = (excludeSelfinteractionQ ? i + 1 : 0); j < mean_.cols(); ++j) {
                        out << Key << j << Value
                            << Flow << BeginSeq
                            << mean()(i, j)
                            << standardError()(i, j)
                            << EndSeq;
                    }
                    out << EndMap;
                } else { // ColumnVector
                    out << Flow << BeginSeq
                        << mean()(i)
                        << standardError()(i)
                        << EndSeq;
                }
            }
            out << EndMap;
        }

        static bool getYAMLMatrixDimensions(const YAML::Node &node, size_t &rows, size_t &cols) {
            if (node.Type() == YAML::NodeType::Map &&
                node.begin()->second.Type() == YAML::NodeType::Map &&
                node.begin()->second.begin()->second.Type() == YAML::NodeType::Sequence) { // Matrix
                rows = node.size();
                cols = node.begin()->second.size();
                return true;
            } else if (node.Type() == YAML::NodeType::Map &&
                       node.begin()->second.Type() == YAML::NodeType::Sequence) { // Vector
                rows = rows = node.size();
                cols = 1;
                return true;
            } else
                return false;
        }

        static MeanErrorPair<Derived> fromYaml(const YAML::Node &node, bool excludeSelfinteractionQ = false) {
            Derived mean, error;
            size_t rows, cols;

            if (!getYAMLMatrixDimensions(node, rows, cols)) throw std::runtime_error("Wrong type for MeanErrorPair.");

            if (excludeSelfinteractionQ) {

                mean.resize(rows + 1, cols + 1);
                error.resize(rows + 1, cols + 1);
                mean.setZero();
                error.setZero();

                for (int i = 0; i < rows; ++i) {
                    for (int j = i + 1; j < cols + 1; ++j) {
                        assert(j < cols + 1 && "To many elements in row.");
                        mean(i, j) = node[i][j][0].as<double>();
                        error(i, j) = node[i][j][1].as<double>();
                    }
                }

            } else {
//TODO also print with self interactions as upper triangular matrix?

                if (cols > 1) {
                    mean.resize(rows, cols);
                    error.resize(rows, cols);
                    mean.setZero();
                    error.setZero();

                    for (int i = 0; i < rows; ++i) {
                        for (int j = 0; j < node[i].size(); ++j) {
                            assert(j < cols && "To many elements in row.");
                            mean(i, j) = node[i][j][0].as<double>();
                            error(i, j) = node[i][j][1].as<double>();
                        }
                    }
                } else {
                    mean.resize(rows, 1);
                    error.resize(rows, 1);
                    mean.setZero();
                    error.setZero();
                    for (int i = 0; i < rows; ++i) {
                        mean(i, 0) = node[i][0].as<double>();
                        error(i, 0) = node[i][1].as<double>();
                    }
                }
            }

            return MeanErrorPair<Derived>(mean, error);
        };

    private:
        void initialize(const Derived &sample) {
            initializedQ_ = true;
            squaredTotalWeightInitializedQ_ = true;
            mean_ = sample;
            lastMean_ = mean_;
            unnormalisedVariance_ = Derived::Zero(sample.rows(), sample.cols());
            lastUnnormalisedVariance_ = unnormalisedVariance_;

            cwiseMin_ = sample;
            cwiseMax_ = sample;
        }

        //https://en.wikipedia.org/wiki/Variance#Population_variance_and_sample_variance
        Derived unbiasedSampleVariance() const { // includes Bessel's correction
            if (totalWeight_ <= 1)
                return Derived::Zero(unnormalisedVariance_.rows(), unnormalisedVariance_.cols());
            else
                return (unnormalisedVariance_ / (totalWeight_ - 1));
        }

        Derived biasedSampleVariance() const { // population variance
            if (totalWeight_ <= 0)
                return Derived::Zero(unnormalisedVariance_.rows(), unnormalisedVariance_.cols());
            else
                return (unnormalisedVariance_ / totalWeight_);
        }

        Derived sampleReliabilityVariance() const {
            if(!std::is_integral<WeightType>::value)
                std::cerr << "The squared total weight was not initialized properly and might be unreliable." << std::endl;
            
            if (totalWeight_ <= 1)
                return Derived::Zero(unnormalisedVariance_.rows(), unnormalisedVariance_.cols());
            else
                return unnormalisedVariance_ / double(totalWeight_ - double(squaredTotalWeight_) / totalWeight_);
        }

        Derived biasedStandardDeviation() const { return biasedSampleVariance().array().sqrt(); }

        Derived unbiasedSampleStandardDeviation() const { return unbiasedSampleVariance().array().sqrt(); }

        Derived calculateVarianceFromStandardError(const Derived &error, WeightType N) const {
            return N * (N - 1) * error.array().pow(2);//TODO CHECK!!
        }

        bool initializedQ_, squaredTotalWeightInitializedQ_;
        WeightType totalWeight_, squaredTotalWeight_;
        Derived mean_, lastMean_, unnormalisedVariance_, lastUnnormalisedVariance_, cwiseMin_, cwiseMax_;
    };
}

#endif //INPSIGHTS_STATISTICS_H
