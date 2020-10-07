// Edited by Leonard Reuter on 26.06.19.
// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <BestMatch.h>
#include <SOAPSettings.h>

template <>
bool BestMatch::AscendingMetricResult::operator<(const BestMatch::Result<true> &rhs) {
    if(metric == rhs.metric)
        for (Eigen::Index i = 0; i < permutation.indices().size(); ++i) {
            if(permutation.indices()[i] != rhs.permutation.indices()[i])
                return permutation.indices()[i] < rhs.permutation.indices()[i];
        }

    return metric < rhs.metric;
}

template <>
bool BestMatch::DescendingMetricResult::operator<(const BestMatch::Result<false> &rhs) {
    auto numericalPrecisionEpsilon = SOAP::General::settings.comparisonEpsilon.get();

    if(abs(metric - rhs.metric) < numericalPrecisionEpsilon)
        for (Eigen::Index i = 0; i < permutation.indices().size(); ++i) {
            if(permutation.indices()[i] != rhs.permutation.indices()[i])
                return permutation.indices()[i] < rhs.permutation.indices()[i];
        }

    return metric > rhs.metric;
}
