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

TEST_F(AStatisticsTest, Constructor) {
    unsigned N = 3;
    auto mean = mat2;
    auto standardError = mat1 * 1/std::sqrt(N);
    auto cwiseMin = mat1;
    auto cwiseMax = mat3;
    auto stats = Statistics::RunningStatistics<Eigen::MatrixXd>(mean, standardError, cwiseMin, cwiseMax, N, N*N);

    stats.add(mat2);

    auto expectedStats = Statistics::RunningStatistics<Eigen::MatrixXd>();
    expectedStats.add(mat1);
    expectedStats.add(mat2);
    expectedStats.add(mat3);
    expectedStats.add(mat2);

    ASSERT_TRUE(stats.mean().isApprox(expectedStats.mean()));
    ASSERT_TRUE(stats.standardDeviation().isApprox(expectedStats.standardDeviation()));
    ASSERT_TRUE(stats.cwiseMin().isApprox(expectedStats.cwiseMin()));
    ASSERT_TRUE(stats.cwiseMax().isApprox(expectedStats.cwiseMax()));
    ASSERT_EQ(stats.getTotalWeight(), expectedStats.getTotalWeight());
    ASSERT_EQ(stats.getSquaredTotalWeight(), expectedStats.getSquaredTotalWeight());
}

TEST_F(AStatisticsTest, Matrix){
    Statistics::RunningStatistics<Eigen::MatrixXd> stats;
    stats.add(mat1);
    stats.add(mat2);
    stats.add(mat3);

    ASSERT_TRUE(stats.mean().isApprox(mat2));
    ASSERT_TRUE(stats.standardDeviation().isApprox(mat1));
    ASSERT_TRUE(stats.cwiseMin().isApprox(mat1));
    ASSERT_TRUE(stats.cwiseMax().isApprox(mat3));
}

TEST_F(AStatisticsTest, Vector){
    Statistics::RunningStatistics<Eigen::VectorXd> stats;
    stats.add(vec1);
    stats.add(vec2);
    stats.add(vec3);

    ASSERT_TRUE(stats.mean().isApprox(vec2));
    ASSERT_TRUE(stats.standardDeviation().isApprox(vec1));
    ASSERT_TRUE(stats.cwiseMin().isApprox(vec1));
    ASSERT_TRUE(stats.cwiseMax().isApprox(vec3));
}

TEST_F(AStatisticsTest, UnsignedIntegerWeights){
    // unsigned is the default template parameter
    Statistics::RunningStatistics<Eigen::VectorXd> stats;

    Eigen::VectorXd vec1(2);
    vec1 << -2,-4;
    Eigen::VectorXd vec2(2);
    vec2 << 1,2;


    stats.add(vec1,1);
    stats.add(vec2,2);

    Eigen::VectorXd expectedMean(2);
    expectedMean << 0,0;
    Eigen::VectorXd expectedVariance(2);
    expectedVariance << 3,12;

    ASSERT_TRUE(stats.mean().isApprox(expectedMean));
    ASSERT_TRUE(stats.variance().isApprox(expectedVariance));
    ASSERT_TRUE(stats.cwiseMin().isApprox(vec1));
    ASSERT_TRUE(stats.cwiseMax().isApprox(vec2));
    
    // Reset
    stats.reset();

    stats.add(vec1,1);
    stats.add(vec2,1);
    stats.add(vec2,1);

    ASSERT_TRUE(stats.mean().isApprox(expectedMean));
    ASSERT_TRUE(stats.variance().isApprox(expectedVariance));
    ASSERT_TRUE(stats.cwiseMin().isApprox(vec1));
    ASSERT_TRUE(stats.cwiseMax().isApprox(vec2));
}

TEST_F(AStatisticsTest, DoublePrecisionWeights){
    Statistics::RunningStatistics<Eigen::VectorXd,double> stats;

    Eigen::VectorXd vec1(2);
    vec1 << -2,-4;
    Eigen::VectorXd vec2(2);
    vec2 << 1,2;

    stats.add(vec1,1.0);
    stats.add(vec2,2.0);

    Eigen::VectorXd expectedMean(2);
    expectedMean << 0,0;
    Eigen::VectorXd expectedVariance(2);
    expectedVariance << 3,12;

    ASSERT_TRUE(stats.mean().isApprox(expectedMean));
    ASSERT_TRUE(stats.variance().isApprox(expectedVariance));
    ASSERT_TRUE(stats.cwiseMin().isApprox(vec1));
    ASSERT_TRUE(stats.cwiseMax().isApprox(vec2));

    // Reset
    stats.reset();

    stats.add(vec1);
    stats.add(vec2);
    stats.add(vec2);

    ASSERT_TRUE(stats.mean().isApprox(expectedMean));
    ASSERT_TRUE(stats.variance().isApprox(expectedVariance));
    ASSERT_TRUE(stats.cwiseMin().isApprox(vec1));
    ASSERT_TRUE(stats.cwiseMax().isApprox(vec2));
}

TEST_F(AStatisticsTest, MinMax){
    Statistics::RunningStatistics<Eigen::VectorXd> stats;
    
    Eigen::VectorXd v1(2),v2(2);
    v1 << -1, 2;
    v2 << 1,-2;
    
    stats.add(v1);
    stats.add(v2);

    Eigen::VectorXd min(2),max(2);
    min << -1,-2;
    max <<  1, 2;
    
    ASSERT_TRUE(stats.cwiseMin().isApprox(min));
    ASSERT_TRUE(stats.cwiseMax().isApprox(max));
    
}

TEST_F(AStatisticsTest, SingleValue){
    Statistics::RunningStatistics<Eigen::Matrix<double,1,1>> stats;

    Eigen::VectorXd v1(1),v2(1);
    v1 << -1;
    v2 << 1;

    stats.add(v1);
    stats.add(v2);

    Eigen::VectorXd min(1),max(1);
    min << -1;
    max <<  1;

    ASSERT_TRUE(stats.cwiseMin().isApprox(min));
    ASSERT_TRUE(stats.cwiseMax().isApprox(max));

}

TEST_F(AStatisticsTest, YAMLConversion){
    Statistics::RunningStatistics<Eigen::MatrixXd> stats;
    stats.add(mat1);
    stats.add(mat2);
    stats.add(mat3);

    YAML::Emitter out;
    out << stats;

    auto loaded = YAML::Load(out.c_str());
    auto readStats = loaded.as<Statistics::RunningStatistics<Eigen::MatrixXd>>();

    ASSERT_TRUE(readStats.mean().isApprox(stats.mean()));
    ASSERT_TRUE(readStats.standardDeviation().isApprox(stats.standardDeviation()));
    ASSERT_TRUE(readStats.cwiseMin().isApprox(stats.cwiseMin()));
    ASSERT_TRUE(readStats.cwiseMax().isApprox(stats.cwiseMax()));
    ASSERT_EQ(readStats.getTotalWeight(), stats.getTotalWeight());
    ASSERT_EQ(readStats.getSquaredTotalWeight(), stats.getSquaredTotalWeight());
}
