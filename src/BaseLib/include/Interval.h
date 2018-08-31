//
// Created by Michael Heuer on 29.08.18.
//

#ifndef AMOLQCPP_INTERVAL_H
#define AMOLQCPP_INTERVAL_H

#include <cassert>

class Interval {
public:
    Interval() = default;

    Interval(long start, long n)
            : start_(start),n_(n) {
        assert(start_ >= 0);
        assert(n_ >= 0);
    }

    bool operator==(const Interval& other)const{
        return (start_ == other.start_) && (n_ == other.n_);
    }

    explicit Interval(long idx)
            : Interval({idx,1}){}

    long start() const { return start_;};
    long end() const { return start()+n_;};
    long numberOfEntities() const { return n_;};

private:
    long start_,n_;
};

#endif //AMOLQCPP_INTERVAL_H
