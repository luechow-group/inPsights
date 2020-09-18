// Edited by Leonard Reuter on 26.06.19.
// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_BESTMATCH_H
#define INPSIGHTS_BESTMATCH_H

#include <Eigen/Core>

namespace BestMatch {
    template<bool ascending = true>
    struct Result {
        double metric;
        Eigen::PermutationMatrix<Eigen::Dynamic> permutation;

        bool operator<(const Result &rhs);
    };

    using AscendingMetricResult = Result<true>;
    using DescendingMetricResult = Result<false>;
};

#endif //INPSIGHTS_BESTMATCH_H
