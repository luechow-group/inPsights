//
// Created by Michael Heuer on 14.09.18.
//

#ifndef AMOLQCPP_STATISTICS_H
#define AMOLQCPP_STATISTICS_H

#include <Eigen/Core>
#include <yaml-cpp/yaml.h>

namespace Statistics {

    template<typename Derived, typename WeightType = unsigned>
    class RunningStatistics {

    public:
        RunningStatistics()
                : initializedQ_(false),totalWeight_(0),squaredTotalWeight_(0),
                mean_(),lastMean_(),
                unnormalisedVariance_(),lastUnnormalisedVariance_() {}

        void reset() {
            initializedQ_ = false;
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
            squaredTotalWeight_ += w*w;

            if (!initializedQ_) {
                initialize(sample);
            } else {
                Derived delta = sample - lastMean_;
                mean_ = lastMean_ + (delta * w/totalWeight_).matrix();
                unnormalisedVariance_ = lastUnnormalisedVariance_ + w*(delta.array() * (sample - mean_).array()).matrix();

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
            return  unbiasedSampleVariance();
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
            for (Eigen::Index i = 0; i < mean_.rows()- (excludeSelfinteractionQ? 1 : 0); ++i) {
                out << Key << i << Value;

                if (mean_.cols() > 1) { // Matrix
                    out << BeginMap;
                    for (Eigen::Index j = (excludeSelfinteractionQ ? i+1 : 0); j < mean_.cols(); ++j) {
                        out << Key << j << Value
                            << Flow << BeginSeq
                            << mean()(i, j);
                        if(getTotalWeight() > 1) out << standardError()(i, j);
                        out << EndSeq;
                    }
                    out << EndMap;
                } else { // ColumnVector
                    out << Flow << BeginSeq
                        << mean()(i);
                    if(getTotalWeight() > 1) out << standardError()(i);
                    out << EndSeq;
                }
            }
            out << EndMap;
        }


    private:
        void initialize(const Derived &sample){
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
            if(totalWeight_ <= 1)
                return Derived::Zero(unnormalisedVariance_.rows(),unnormalisedVariance_.cols());
            else
                return (unnormalisedVariance_ / (totalWeight_ - 1));
        }

        Derived biasedSampleVariance() const { // population variance
            if(totalWeight_ <= 0)
                return Derived::Zero(unnormalisedVariance_.rows(),unnormalisedVariance_.cols());
            else
                return (unnormalisedVariance_ / totalWeight_);
        }

        Derived sampleReliabilityVariance() const {
            if(totalWeight_ <= 1)
                return Derived::Zero(unnormalisedVariance_.rows(),unnormalisedVariance_.cols());
            else
                return unnormalisedVariance_ / double(totalWeight_ - double(squaredTotalWeight_)/totalWeight_);
        }

        Derived biasedStandardDeviation() const { return biasedSampleVariance().array().sqrt(); }

        Derived unbiasedSampleStandardDeviation() const { return unbiasedSampleVariance().array().sqrt(); }

        bool initializedQ_;
        WeightType totalWeight_,squaredTotalWeight_;
        Derived mean_, lastMean_, unnormalisedVariance_, lastUnnormalisedVariance_, cwiseMin_, cwiseMax_;
    };
}

#endif //AMOLQCPP_STATISTICS_H
