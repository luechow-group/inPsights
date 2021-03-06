// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_CLUSTERINGMETRIC_H
#define INPSIGHTS_CLUSTERINGMETRIC_H


#include <EnvironmentBasedMetric.h>
#include <DistanceBasedMetric.h>
#include "Cluster.h"
#include "ErrorHandling.h"
#include "Maximum.h"

enum MetricType {Environmental, Spatial};
enum SelectionType {Global, Local};

// templat für environment oder distance based
namespace ClusteringMetric{

    //BestMatch::AscendingMetricResult allResults(const Cluster &g1, const Cluster &g2); //optional

    template <MetricType metric, SelectionType selection>
    BestMatch::AscendingMetricResult bestMatchResult(const Cluster &g1, const Cluster &g2);

    // Forward declaration of the template
    template<> BestMatch::AscendingMetricResult
    bestMatchResult<Spatial, Global> (const Cluster &g1, const Cluster &g2);
    template<> BestMatch::AscendingMetricResult
    bestMatchResult<Spatial, Local> (const Cluster &g1, const Cluster &g2);

    template<> BestMatch::AscendingMetricResult
    bestMatchResult<Environmental, Global> (const Cluster &g1, const Cluster &g2);
    template<> BestMatch::AscendingMetricResult
    bestMatchResult<Environmental, Local> (const Cluster &g1, const Cluster &g2);


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
