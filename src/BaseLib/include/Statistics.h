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
                : initializedQ_(false),wSum_(0),wSum2_(0),
                X_(),Xold_(),
                M2_(),M2old_() {}

        void reset() {
            initializedQ_ = false;
            wSum_ = 0;
            wSum2_ = 0;
        }

        void initialize(const Derived &sample){
            initializedQ_ = true;
            X_ = sample;
            Xold_ = X_;
            M2old_ = Derived::Zero(sample.rows(), sample.cols());
        }
        
        void add(const Derived &sample, WeightType w = 1) {
            wSum_ += w;
            wSum2_ += w*w;

            if (!initializedQ_) {
                initialize(sample);
            } else {
                Derived delta = sample - Xold_;
                X_ = Xold_ + (delta * w/wSum_).matrix();
                M2_ = M2old_ + w*(delta.array() * (sample - X_).array()).matrix();

                Xold_ = X_;
                M2old_ = M2_;
            }
        }

        Derived mean() {
            assert(wSum_ >= 1 && "The number of samples must be at least one.");
            return X_.matrix();
        }

        Derived variance(){
            return  unbiasedSampleVariance();
        }

        Derived standardDeviation(){
            return unbiasedSampleStandardDeviation();
        }

    private:
        //https://en.wikipedia.org/wiki/Variance#Population_variance_and_sample_variance
        Derived unbiasedSampleVariance() { // includes Bessel's correction
            assert(wSum_ >= 2 && "The number of samples must be at least two.");
            return (M2_ / (wSum_ - 1));
        }

        Derived biasedSampleVariance() { // population variance
            assert(wSum_ >= 2 && "The number of samples must be at least two.");
            return (M2_ / wSum_);
        }

        Derived sampleReliabilityVariance() {
            assert(wSum_ >= 2 && "The number of samples must be at least two.");
            return M2_ / double(wSum_ - double(wSum2_)/wSum_);
        }

        Derived biasedStandardDeviation() { return biasedSampleVariance().array().sqrt(); }

        Derived unbiasedSampleStandardDeviation() { return unbiasedSampleVariance().array().sqrt(); }



        bool initializedQ_;
        WeightType wSum_,wSum2_;
        Derived X_, Xold_, M2_, M2old_;
    };
}


//sample_reliability_variance = S / (wSum - wSum2/wSum)

#endif //AMOLQCPP_STATISTICS_H
