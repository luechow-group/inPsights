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

#ifndef INPSIGHTS_GENERALSTATISTICS_H
#define INPSIGHTS_GENERALSTATISTICS_H

#include "Reference.h"
#include "Sample.h"
#include "Group.h"
#include <Statistics.h>

namespace GeneralStatistics{
    struct Result { //TODO make class
        Result() = default;

        //TODO better names
        Result(SingleValueStatistics valueStats,
                SingleValueStatistics EelStats,
                SingleValueStatistics TeStats,
                SingleValueStatistics VeeStats,
                SingleValueStatistics VenStats,
                double Vnn)
                :
                valueStats_(std::move(valueStats)),
                EStats_(std::move(EelStats)),
                TeStats_(std::move(TeStats)),
                VeeStats_(std::move(VeeStats)),
                VenStats_(std::move(VenStats)),
                Vnn_(Vnn) {};

        SingleValueStatistics valueStats_, EStats_, TeStats_,VeeStats_, VenStats_;
        double Vnn_;
    };

    Result calculate(Group &maxima,
                     std::vector<Sample> &samples,
                     const AtomsVector& atoms);
}

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<GeneralStatistics::Result> {
        static Node encode(const GeneralStatistics::Result &rhs);
        static bool decode(const Node &node, GeneralStatistics::Result &rhs);
    };
    Emitter &operator<<(Emitter &out, const GeneralStatistics::Result &p) ;
}

#endif //INPSIGHTS_GENERALSTATISTICS_H
