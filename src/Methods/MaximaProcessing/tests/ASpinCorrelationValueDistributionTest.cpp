/* Copyright 2020 Michael Heuer
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
