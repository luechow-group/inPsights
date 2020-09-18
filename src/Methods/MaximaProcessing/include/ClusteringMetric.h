// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CLUSTERINGMETRIC_H
#define INPSIGHTS_CLUSTERINGMETRIC_H


#include <EnvironmentBasedMetric.h>
#include <DistanceBasedMetric.h>
#include "Cluster.h"
#include "ErrorHandling.h"
#include "Reference.h"

enum MetricType {Environmental, Spatial};
enum SelectionType {Global, Local};

// templat für environment oder distance based
namespace ClusteringMetric{

    //BestMatch::AscendingMetricResult allResults(const Cluster &g1, const Cluster &g2); //optional

    template <MetricType metric, SelectionType selection>
    BestMatch::AscendingMetricResult bestMatchResult(const Cluster &g1, const Cluster &g2){
        throw NotImplemented();
        return {0, Eigen::PermutationMatrix<Eigen::Dynamic>(0)};
    }

    template<> BestMatch::AscendingMetricResult bestMatchResult<Spatial, Global> (const Cluster &g1, const Cluster &g2) {
        return Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
                g1.representative()->maximum().positionsVector(),
                g2.representative()->maximum().positionsVector());
    }

    template<> BestMatch::AscendingMetricResult bestMatchResult<Spatial, Local> (const Cluster &g1, const Cluster &g2) {
        auto g1ElectronsCount = g1.getSelectedElectronsCount();
        auto g2ElectronsCount = g2.getSelectedElectronsCount();

        if (g1ElectronsCount == g2ElectronsCount) {
            return Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
                    g1.representative()->maximum().head(g1ElectronsCount).positionsVector(),
                    g2.representative()->maximum().head(g2ElectronsCount).positionsVector());
        }
        //TODO Should we throw an exception instead of return max here?
        return {std::numeric_limits<double>::max(), Eigen::PermutationMatrix<Eigen::Dynamic>(0)};
    }

    template <MetricType metric, SelectionType selection>
    double distance(const Cluster &g1, const Cluster &g2){
        return bestMatchResult<metric, selection>(g1, g2).metric;
    }

    template <MetricType metric, SelectionType selection>
    bool similarQ(const Cluster &g1, const Cluster &g2, double threshold){
        return distance<metric, selection>(g1, g2) < threshold;
    }
};



#endif //INPSIGHTS_CLUSTERINGMETRIC_H
