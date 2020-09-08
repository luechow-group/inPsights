// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <SpinCorrelationValueHistogram.h>

SpinCorrelationValueHistogram::SpinCorrelationValueHistogram(Eigen::Index oneSidedNonzeroBinCount)
:
oneSidedNonzeroBinCount_(oneSidedNonzeroBinCount),
binCount_(oneSidedNonzeroBinCount_ * 2 + 1),
bins_(Eigen::VectorXd::Zero(binCount_)){};

Eigen::Index SpinCorrelationValueHistogram::calculateBinIndex(double spinCorrelation){
    auto binLength = 2.0 / static_cast<double>(binCount_);
    Eigen::Index binIndex = std::ceil(std::abs(spinCorrelation) / binLength - 0.5);
    return spinCorrelation >= 0 ? oneSidedNonzeroBinCount_ + binIndex : oneSidedNonzeroBinCount_ - binIndex;
};

void SpinCorrelationValueHistogram::addSpinStatistic(const TriangularMatrixStatistics& spinCorrelations) {
    assert(spinCorrelations.mean().minCoeff() >= -1.0);
    assert(spinCorrelations.mean().maxCoeff() <= 1.0);

    for (Eigen::Index i = 0; i < spinCorrelations.rows()-1; ++i) {
        for (Eigen::Index j = i+1; j < spinCorrelations.cols(); ++j) {
            auto spinCorrelationValue = spinCorrelations.mean()(i,j);
            bins_[calculateBinIndex(spinCorrelationValue)]  += static_cast<double>(spinCorrelations.getTotalWeight());
        }
    }
}

Eigen::VectorXd SpinCorrelationValueHistogram::getHistogramVector(){
    return bins_;
};
