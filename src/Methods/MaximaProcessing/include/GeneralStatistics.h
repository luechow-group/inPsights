// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_GENERALSTATISTICS_H
#define INPSIGHTS_GENERALSTATISTICS_H

#include "Maximum.h"
#include "Sample.h"
#include "Cluster.h"
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

    Result calculate(Cluster &maxima,
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
