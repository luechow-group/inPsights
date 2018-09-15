//
// Created by Michael Heuer on 15.09.18.
//

#include <gtest/gtest.h>
#include <Statistics.h>
#include <iostream>

using namespace testing;
using namespace Eigen;


class AStatisticsTest : public Test {
public:
    Eigen::MatrixXd mat1, mat2, mat3;
    Eigen::VectorXd vec1, vec2, vec3;

    void SetUp() override {
        mat1.resize(2,2);
        mat1 << 2,4,6,8;

        mat2.resize(2,2);
        mat2 << 4,8,12,16;

        mat3.resize(2,2);
        mat3 << 6,12,18,24;

        vec1.resize(4);
        vec1 << 2,4,6,8;

        vec2.resize(4);
        vec2 << 4,8,12,16;

        vec3.resize(4);
        vec3 << 6,12,18,24;
    }
};

TEST_F(AStatisticsTest, Matrix){
    Statistics::RunningStatistics<Eigen::MatrixXd> meanAndVariance;
    meanAndVariance.add(mat1);
    meanAndVariance.add(mat2);
    meanAndVariance.add(mat3);

    ASSERT_TRUE(meanAndVariance.mean().isApprox(mat2));
    ASSERT_TRUE(meanAndVariance.standardDeviation().isApprox(mat1));
}

TEST_F(AStatisticsTest, Vector){
    Statistics::RunningStatistics<Eigen::VectorXd> meanAndVariance;
    meanAndVariance.add(vec1);
    meanAndVariance.add(vec2);
    meanAndVariance.add(vec3);

    ASSERT_TRUE(meanAndVariance.mean().isApprox(vec2));
    ASSERT_TRUE(meanAndVariance.standardDeviation().isApprox(vec1));
}

TEST_F(AStatisticsTest, UnsignedIntegerWeights){
    // unsigned is the default template parameter
    Statistics::RunningStatistics<Eigen::VectorXd> meanAndVariance;

    Eigen::VectorXd vec1(2);
    vec1 << -2,-4;
    Eigen::VectorXd vec2(2);
    vec2 << 1,2;


    meanAndVariance.add(vec1,1);
    meanAndVariance.add(vec2,2);

    Eigen::VectorXd expectedMean(2);
    expectedMean << 0,0;
    Eigen::VectorXd expectedVariance(2);
    expectedVariance << 3,12;

    ASSERT_TRUE(meanAndVariance.mean().isApprox(expectedMean));
    ASSERT_TRUE(meanAndVariance.variance().isApprox(expectedVariance));

    // Reset
    meanAndVariance.reset();

    meanAndVariance.add(vec1,1);
    meanAndVariance.add(vec2,1);
    meanAndVariance.add(vec2,1);

    ASSERT_TRUE(meanAndVariance.mean().isApprox(expectedMean));
    ASSERT_TRUE(meanAndVariance.variance().isApprox(expectedVariance));
}

TEST_F(AStatisticsTest, DoublePrecisionWeights){
    Statistics::RunningStatistics<Eigen::VectorXd,double> meanAndVariance;

    Eigen::VectorXd vec1(2);
    vec1 << -2,-4;
    Eigen::VectorXd vec2(2);
    vec2 << 1,2;

    meanAndVariance.add(vec1,1.0);
    meanAndVariance.add(vec2,2.0);

    Eigen::VectorXd expectedMean(2);
    expectedMean << 0,0;
    Eigen::VectorXd expectedVariance(2);
    expectedVariance << 3,12;

    ASSERT_TRUE(meanAndVariance.mean().isApprox(expectedMean));
    ASSERT_TRUE(meanAndVariance.variance().isApprox(expectedVariance));

    // Reset
    meanAndVariance.reset();

    meanAndVariance.add(vec1);
    meanAndVariance.add(vec2);
    meanAndVariance.add(vec2);

    ASSERT_TRUE(meanAndVariance.mean().isApprox(expectedMean));
    ASSERT_TRUE(meanAndVariance.variance().isApprox(expectedVariance));
}