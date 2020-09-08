// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <spdlog/spdlog.h>
#include <SpinCorrelationValueHistogram.h>

using namespace testing;

class ASpinCorrelationValueDistributionTest : public ::testing::Test {
public:

    TriangularMatrixStatistics spinCorrelationsStatistic;
    void SetUp() override {
        Eigen::MatrixXd spinCorrelations;
        spinCorrelations.resize(6,6);
        spinCorrelations <<
        0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        0.0, 0.0, 0.5, 0.5, 0.5, 0.5,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,-0.5,-0.5,
        0.0, 0.0, 0.0, 0.0, 0.0,-1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

        spinCorrelationsStatistic.add(spinCorrelations);
    }
};

TEST_F(ASpinCorrelationValueDistributionTest, fiveBins) {
    SpinCorrelationValueHistogram distribution(2);
    distribution.addSpinStatistic(spinCorrelationsStatistic);

    auto hist = distribution.getHistogramVector();
    Eigen::VectorXd expected(5);
    expected << 1,2,3,4,5;
    ASSERT_EQ(hist, expected);

    // add twice
    distribution.addSpinStatistic(spinCorrelationsStatistic);
    hist = distribution.getHistogramVector();
    ASSERT_EQ(hist, 2*expected);
}

TEST_F(ASpinCorrelationValueDistributionTest, calculateBinIndex) {
    SpinCorrelationValueHistogram distribution(2);
    distribution.addSpinStatistic(spinCorrelationsStatistic);

    ASSERT_EQ(distribution.calculateBinIndex(1.0), 4);
    ASSERT_EQ(distribution.calculateBinIndex(0.8), 4);
    ASSERT_EQ(distribution.calculateBinIndex(0.6), 3);
    ASSERT_EQ(distribution.calculateBinIndex(0.4), 3);
    ASSERT_EQ(distribution.calculateBinIndex(0.2), 2);
    ASSERT_EQ(distribution.calculateBinIndex(0.0), 2);
    ASSERT_EQ(distribution.calculateBinIndex(-0.2), 2);
    ASSERT_EQ(distribution.calculateBinIndex(-0.4), 1);
    ASSERT_EQ(distribution.calculateBinIndex(-0.6), 1);
    ASSERT_EQ(distribution.calculateBinIndex(-0.8), 0);
    ASSERT_EQ(distribution.calculateBinIndex(-1.0), 0);
}
