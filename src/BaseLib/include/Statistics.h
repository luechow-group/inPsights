//
// Created by Michael Heuer on 14.09.18.
//

#ifndef AMOLQCPP_STATISTICS_H
#define AMOLQCPP_STATISTICS_H

#include <Eigen/Core>

namespace Statistics {

    template<typename Derived>
    class RunningStatistics {

    public:
        RunningStatistics()
                : count_(0), X_() {}

        explicit RunningStatistics(const Derived &firstSample)
                :
                count_(1),
                X_(firstSample),
                Xold_(X_),
                M2_(),
                M2old_(Derived::Zero(firstSample.rows(), firstSample.cols())) {};

        void reset() {
            count_ = 0;
        }

        void add(const Derived &sample) {
            count_++;
            if (count_ == 1) {
                X_ = sample;
                Xold_ = X_;
                M2old_ = Derived::Zero(sample.rows(), sample.cols());
            } else {
                Derived delta = sample - Xold_;
                X_ = Xold_ + (delta / count_).matrix();
                M2_ = M2old_ + (delta.array() * (sample - X_).array()).matrix();

                Xold_ = X_;
                M2old_ = M2_;
            }
        }

        Derived mean() {
            assert(count_ >= 1 && "The number of samples must be at least one.");
            return X_.matrix();
        }

        Derived variance() {
            assert(count_ >= 2 && "The number of samples must be at least two.");
            return (M2_ / (count_ - 1));
        }

        Derived sampleVariance() {
            assert(count_ >= 2 && "The number of samples must be at least two.");
            return (M2_ / count_);
        }

        Derived standardDeviation() { return variance().array().sqrt(); }

        Derived sampleStandardDeviation() { return sampleVariance().array().sqrt(); }


    private:
        long count_;
        Derived X_, Xold_, M2_, M2old_;
    };
}

#endif //AMOLQCPP_STATISTICS_H
