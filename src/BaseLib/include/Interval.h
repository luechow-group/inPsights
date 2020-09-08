// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_INTERVAL_H
#define INPSIGHTS_INTERVAL_H

#include <cassert>

class Interval {
public:
    Interval() = default;

    Interval(long start, long n)
            : start_(start),n_(n) {
        assert(start_ >= 0 && "Start index must be greater non-negative.");
        assert(n_ >= 0 && "Interval length must be non-negative");
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

#endif //INPSIGHTS_INTERVAL_H
