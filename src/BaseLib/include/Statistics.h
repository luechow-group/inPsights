//
// Created by Michael Heuer on 14.09.18.
//

#ifndef AMOLQCPP_STATISTICS_H
#define AMOLQCPP_STATISTICS_H

#include <Eigen/Core>

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
        }

        void initialize(const Derived &sample){
            initializedQ_ = true;
            mean_ = sample;
            lastMean_ = mean_;
            unnormalisedVariance_ = Derived::Zero(sample.rows(), sample.cols());
            lastUnnormalisedVariance_ = unnormalisedVariance_;
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
            }
        }

        Derived mean() {
            assert(totalWeight_ >= 1 && "The number of samples must be at least one.");
            return mean_.matrix();
        }

        Derived variance(){
            return  unbiasedSampleVariance();
        }

        Derived standardDeviation(){
            return unbiasedSampleStandardDeviation();
        }

        WeightType getTotalWeight() const {
            return totalWeight_;
        }

    private:
        //https://en.wikipedia.org/wiki/Variance#Population_variance_and_sample_variance
        Derived unbiasedSampleVariance() { // includes Bessel's correction
            assert(totalWeight_ >= 2 && "The number of samples must be at least two.");
            return (unnormalisedVariance_ / (totalWeight_ - 1));
        }

        Derived biasedSampleVariance() { // population variance
            assert(totalWeight_ >= 2 && "The number of samples must be at least two.");
            return (unnormalisedVariance_ / totalWeight_);
        }

        Derived sampleReliabilityVariance() {
            assert(totalWeight_ >= 2 && "The number of samples must be at least two.");
            return unnormalisedVariance_ / double(totalWeight_ - double(squaredTotalWeight_)/totalWeight_);
        }

        Derived biasedStandardDeviation() { return biasedSampleVariance().array().sqrt(); }

        Derived unbiasedSampleStandardDeviation() { return unbiasedSampleVariance().array().sqrt(); }

        bool initializedQ_;
        WeightType totalWeight_,squaredTotalWeight_;
        Derived mean_, lastMean_, unnormalisedVariance_, lastUnnormalisedVariance_;
    };
}

#endif //AMOLQCPP_STATISTICS_H
