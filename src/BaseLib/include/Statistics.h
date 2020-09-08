// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_STATISTICS_H
#define INPSIGHTS_STATISTICS_H

#include <Eigen/Core>
#include <yaml-cpp/yaml.h>
#include <iostream>

namespace Statistics {
    template<typename Derived, typename WeightType = unsigned, bool triangularExport = false>
    class RunningStatistics {

    public:
        RunningStatistics()
                :
                initializedQ_(false),
                totalWeight_(0),
                mean_(), lastMean_(),
                unnormalisedVariance_(), lastUnnormalisedVariance_(),
                cwiseMin_(), cwiseMax_() {}

        RunningStatistics(
                Derived mean, Derived standardError,
                Derived cwiseMin, Derived cwiseMax,
                WeightType totalWeight)
                :
                initializedQ_(true),
                totalWeight_(totalWeight),
                mean_(std::move(mean)),
                lastMean_(mean_),
                unnormalisedVariance_(calculateVarianceFromStandardError(standardError, totalWeight)),
                lastUnnormalisedVariance_(unnormalisedVariance_),
                cwiseMin_(cwiseMin),
                cwiseMax_(cwiseMax) {}


        void reset() {
            initializedQ_ = false;
            totalWeight_ = 0;
            mean_.setZero();
            lastMean_.setZero();
            unnormalisedVariance_.setZero();
            lastUnnormalisedVariance_.setZero();
            cwiseMin_.setZero();
            cwiseMax_.setZero();
        }

        void add(const Derived &sample, WeightType w = 1) {
            totalWeight_ += w;
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

        Eigen::Index rows() const { return mean().rows(); };

        Eigen::Index cols() const { return mean().cols(); };

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
        
    private:
        void initialize(const Derived &sample) {
            initializedQ_ = true;
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
                return Derived::Zero(rows(), cols());
            else
                return (unnormalisedVariance_ / (totalWeight_ - 1));
        }

        Derived biasedSampleVariance() const { // population variance
            if (totalWeight_ <= 0)
                return Derived::Zero(rows(), cols());
            else
                return (unnormalisedVariance_ / totalWeight_);
        }

        Derived biasedStandardDeviation() const { return biasedSampleVariance().array().sqrt(); }

        Derived unbiasedSampleStandardDeviation() const { return unbiasedSampleVariance().array().sqrt(); }

        Derived calculateVarianceFromStandardError(const Derived &error, WeightType N) const {
            return N * (N - 1) * error.array().pow(2);
        }

        bool initializedQ_;
        WeightType totalWeight_;
        Derived mean_, lastMean_, unnormalisedVariance_, lastUnnormalisedVariance_, cwiseMin_, cwiseMax_;
    };
}

using SingleValueStatistics = Statistics::RunningStatistics<Eigen::Matrix<double,1,1>, unsigned>;
using VectorStatistics = Statistics::RunningStatistics<Eigen::VectorXd, unsigned>;
using MatrixStatistics = Statistics::RunningStatistics<Eigen::MatrixXd, unsigned>;
using TriangularMatrixStatistics = Statistics::RunningStatistics<Eigen::MatrixXd, unsigned, true>;

namespace YAML {

    template<typename Derived, typename WeightType, bool triangularExport = false>
    void addStatsSequence(
            const Statistics::RunningStatistics<Derived, WeightType, triangularExport> &rhs,
            Node &node,
            Eigen::Index i,
            Eigen::Index j) {
        node[i][j].push_back(rhs.mean()(i, j));
        node[i][j].push_back(rhs.standardError()(i, j));
        node[i][j].push_back(rhs.cwiseMin()(i, j));
        node[i][j].push_back(rhs.cwiseMax()(i, j));
        node[i]["."] = 0; // Trick to get a map
        node[i].remove(".");

    }

    template<typename Derived, typename WeightType, bool triangularExport = false>
    void addStatsSequence(
            const Statistics::RunningStatistics<Derived, WeightType, triangularExport> &rhs,
            Node &node,
            Eigen::Index i) {
        node[i].push_back(rhs.mean()(i));
        node[i].push_back(rhs.standardError()(i));
        node[i].push_back(rhs.cwiseMin()(i));
        node[i].push_back(rhs.cwiseMax()(i));
        node["."] = 0; // Trick to get a map
        node.remove(".");
    }

    template<typename Derived>
    void readStatsSequence(const Node &node,
                           Derived& mean, Derived& standardError, Derived& cwiseMin, Derived& cwiseMax,
                           size_t i, size_t j) {
        typedef typename Derived::Scalar Scalar;
        mean(i, j) = node[i][j][0].as<Scalar>();
        standardError(i, j) = node[i][j][1].as<Scalar>();
        cwiseMin(i, j) = node[i][j][2].as<Scalar>();
        cwiseMax(i, j) = node[i][j][3].as<Scalar>();
    }

    template<typename Derived>
    void readStatsSequence(const Node &node,
                           Derived& mean, Derived& standardError, Derived& cwiseMin, Derived& cwiseMax, size_t i) {
        typedef typename Derived::Scalar Scalar;
        mean(i) = node[i][0].as<Scalar>();
        standardError(i) = node[i][1].as<Scalar>();
        cwiseMin(i) = node[i][2].as<Scalar>();
        cwiseMax(i) = node[i][3].as<Scalar>();
    }


    template<typename Derived, typename WeightType>
    struct convert<Statistics::RunningStatistics<Derived, WeightType, false>> {
        static Node encode(const Statistics::RunningStatistics<Derived, WeightType, false> &rhs) {
            Node node;

            for (Eigen::Index i = 0; i < rhs.rows(); ++i)
                for (Eigen::Index j = 0; j < rhs.cols(); ++j)
                    if(rhs.cols() > 1)
                        addStatsSequence<Derived,WeightType,false>(rhs, node, i, j);
                    else  // column vector
                        addStatsSequence<Derived,WeightType,false>(rhs, node, i);

            node["N"] = rhs.getTotalWeight();
            return node;
        }

        static bool decode(const Node &node, Statistics::RunningStatistics<Derived, WeightType, false> &rhs) {
            if (!node.IsMap()
                && !node.begin()->second.IsMap()
                && !node.begin()->second.begin()->second.IsSequence())
                return false;

            auto rows = node.size() - 1;
            size_t cols;

            if(node.begin()->second.IsMap())
                cols = node.begin()->second.size();
            else if(node.begin()->second.IsSequence())
                cols = 1;
            else
                return false;

            auto N = node["N"].as<WeightType>();

            Derived mean = Eigen::MatrixBase< Derived >::Zero(rows,cols);
            Derived standardError = Eigen::MatrixBase< Derived >::Zero(rows,cols);
            Derived cwiseMin = Eigen::MatrixBase< Derived >::Zero(rows,cols);
            Derived cwiseMax = Eigen::MatrixBase< Derived >::Zero(rows,cols);

            for (size_t i = 0; i < rows; ++i)
                for (size_t j = 0; j < cols; ++j) {
                    if (cols > 1)
                        readStatsSequence<Derived>(node, mean, standardError, cwiseMin, cwiseMax, i, j);
                    else  // column vector
                        readStatsSequence<Derived>(node, mean, standardError, cwiseMin, cwiseMax, i);
                }


            Statistics::RunningStatistics<Derived, WeightType> stats(mean, standardError, cwiseMin, cwiseMax, N);
            rhs = stats;
            return true;
        }
    };

    template<typename Derived, typename WeightType>
    struct convert<Statistics::RunningStatistics<Derived, WeightType, true>> {
        static Node encode(const Statistics::RunningStatistics<Derived, WeightType, true> &rhs) {
            Node node;

            for (Eigen::Index i = 0; i < rhs.rows()-1; ++i)
                for (Eigen::Index j = i+1; j < rhs.cols(); ++j)
                    addStatsSequence<Derived,WeightType,true>(rhs, node, i, j);

            node["N"] = rhs.getTotalWeight();
            return node;
        }

        static bool decode(const Node &node, Statistics::RunningStatistics<Derived, WeightType,true> &rhs) {
            if (!node.IsMap()
                && !node.begin()->second.IsMap()
                && !node.begin()->second.begin()->second.IsSequence())
                return false;

            auto N = node["N"].as<WeightType>();

            auto rows = node.size() - 1; // because of additional "N" key
            auto cols = node.begin()->second.size();

            Derived mean = Derived::Zero(rows+1, cols+1);
            Derived standardError = Derived::Zero(rows+1, cols+1);
            Derived cwiseMin = Derived::Zero(rows+1, cols+1);
            Derived cwiseMax = Derived::Zero(rows+1, cols+1);

            for (size_t i = 0; i < rows; ++i)
                for (size_t j = i + 1; j < cols + 1; ++j) {
                    assert(j < cols + 1 && "To many elements in row.");
                    readStatsSequence<Derived>(node, mean, standardError, cwiseMin, cwiseMax, i, j);
                }

            Statistics::RunningStatistics<Derived, WeightType,true> stats(mean, standardError, cwiseMin, cwiseMax, N);
            rhs = stats;
            return true;
        }
    };

    template<typename Derived, typename WeightType, bool triangularExport = false>
    void emitStats(Emitter &out,
                   const Statistics::RunningStatistics<Derived, WeightType, triangularExport> &rhs,
                   Eigen::Index i, Eigen::Index j) {
        out << Key << j << Value
            << Flow << BeginSeq
            << rhs.mean()(i, j)
            << rhs.standardError()(i, j)
            << rhs.cwiseMin()(i, j)
            << rhs.cwiseMax()(i, j)
            << EndSeq;
    };

    template<typename Derived, typename WeightType, bool triangularExport = false>
    void emitStats(Emitter &out,
                   const Statistics::RunningStatistics<Derived, WeightType, triangularExport> &rhs,
                   Eigen::Index i) {
        out << Key << i << Value
            << Flow << BeginSeq
            << rhs.mean()(i)
            << rhs.standardError()(i)
            << rhs.cwiseMin()(i)
            << rhs.cwiseMax()(i)
            << EndSeq;
    };

    template<typename Derived, typename WeightType>
    Emitter &operator<<(Emitter &out, const Statistics::RunningStatistics<Derived, WeightType,false> &rhs) {

        out << BeginMap;
        for (Eigen::Index i = 0; i < rhs.rows(); ++i) {
            if(rhs.cols() > 1) {
                out << Key << i << Value << BeginMap;
                for (Eigen::Index j = 0; j < rhs.cols(); ++j) {
                    emitStats<Derived, WeightType, false>(out, rhs, i, j);
                }
                out << EndMap;
            } else { // column vector
                emitStats<Derived, WeightType, false>(out, rhs, i);
            }
        }
        out << Key << "N" << Value << rhs.getTotalWeight() << EndMap;
        return out;
    };


    template<typename Derived, typename WeightType>
    Emitter &operator<<(Emitter &out, const Statistics::RunningStatistics<Derived, WeightType,true> &rhs) {

        out << BeginMap;
        for (Eigen::Index i = 0; i < rhs.rows()-1; ++i) {
            out << Key << i << Value << BeginMap;
            for (Eigen::Index j = i+1; j < rhs.cols(); ++j) {
                emitStats<Derived,WeightType,true>(out, rhs, i, j);
            }
            out << EndMap;
        }
        out << Key << "N" << Value << rhs.getTotalWeight() << EndMap;
        return out;
    }


}

#endif //INPSIGHTS_STATISTICS_H
