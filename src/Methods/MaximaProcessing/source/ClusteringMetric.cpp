// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ClusteringMetric.h"

template<> BestMatch::AscendingMetricResult ClusteringMetric::bestMatchResult<Spatial, Global> (const Cluster &g1, const Cluster &g2) {
    return Metrics::Similarity::DistanceBased::compare<Eigen::Infinity, 2>(
            g1.representative()->maximum().positionsVector(),
            g2.representative()->maximum().positionsVector());
}

template<> BestMatch::AscendingMetricResult ClusteringMetric::bestMatchResult<Spatial, Local> (const Cluster &g1, const Cluster &g2) {
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

template<> BestMatch::AscendingMetricResult ClusteringMetric::bestMatchResult<Environmental, Global> (const Cluster &g1, const Cluster &g2) {
    throw NotImplemented();
    auto descedingMetricResult = Metrics::Similarity::EnvironmentBased::compare(
            g1.representative()->spectrum(),
            g2.representative()->spectrum(),
            0, 0, 0);
    return {1.0-descedingMetricResult.metric, descedingMetricResult.permutation};
}

template<> BestMatch::AscendingMetricResult ClusteringMetric::bestMatchResult<Environmental, Local> (const Cluster &g1, const Cluster &g2) {
    throw NotImplemented();
    return {0, Eigen::PermutationMatrix<Eigen::Dynamic>(0)};
}