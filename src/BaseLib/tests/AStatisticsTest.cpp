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

TEST_F(AStatisticsTest, FirstSampleConstructor){
    Statistics::RunningStatistics<Eigen::MatrixXd> meanAndVariance(mat1);
    meanAndVariance.add(mat2);
    meanAndVariance.add(mat3);

    ASSERT_TRUE(meanAndVariance.mean().isApprox(mat2));
    ASSERT_TRUE(meanAndVariance.standardDeviation().isApprox(mat1));
}

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
