// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_SPINCORRELATIONVALUEHISTOGRAM_H
#define INPSIGHTS_SPINCORRELATIONVALUEHISTOGRAM_H

#include <Statistics.h>

class SpinCorrelationValueHistogram {
public:
    SpinCorrelationValueHistogram(Eigen::Index oneSidedNonzeroBinCount = 12);

    Eigen::Index calculateBinIndex(double spinCorrelation);

    void addSpinStatistic(const TriangularMatrixStatistics& spinCorrelations);

    Eigen::VectorXd getHistogramVector();

private:
    Eigen::Index oneSidedNonzeroBinCount_,binCount_;
    Eigen::VectorXd bins_;
};

#endif //INPSIGHTS_SPINCORRELATIONVALUEHISTOGRAM_H
