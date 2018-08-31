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
            : start_(start),end_(start+n) {
        assert(start_ >= 0);
        assert(start_ <= end_);
    }
    explicit Interval(long idx)
            : Interval({idx,1}){}

    bool checkBounds(long numberOfEntities) const {
        return end() <= numberOfEntities;
    }

    long start() const { return start_;};
    long end() const { return end_;};
    long numberOfEntities() const { return end()-start();};

private:
    long start_,end_;
};

#endif //AMOLQCPP_INTERVAL_H
