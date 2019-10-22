/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
