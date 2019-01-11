//
// Created by Michael Heuer on 18.12.18.
//

#ifndef INPSIGHTS_GENERALSTATISTICS_H
#define INPSIGHTS_GENERALSTATISTICS_H

#include "Reference.h"
#include "SimilarReferences.h"
#include "Sample.h"
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
                EelStats_(std::move(EelStats)),
                TeStats_(std::move(TeStats)),
                VeeStats_(std::move(VeeStats)),
                VenStats_(std::move(VenStats)),
                Vnn_(Vnn) {};

        SingleValueStatistics valueStats_, EelStats_, TeStats_,VeeStats_, VenStats_;
        double Vnn_;
    };

    Result calculate(std::vector<Reference> &references,
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
