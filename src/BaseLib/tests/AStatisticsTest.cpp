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
    auto stats = Statistics::RunningStatistics<Eigen::MatrixXd>(mean, standardError, cwiseMin, cwiseMax, N);

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

TEST_F(AStatisticsTest, YAMLMatrixConversion){
    Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,false> stats;
    stats.add(mat1);
    stats.add(mat2);
    stats.add(mat3);

    YAML::Emitter out;
    out << stats;

    auto expectedString = "0:\n  0: [4, 1.154700538379252, 2, 6]\n  1: [8, 2.309401076758503, 4, 12]\n"\
                          "1:\n  0: [12, 3.464101615137755, 6, 18]\n  1: [16, 4.618802153517007, 8, 24]\n"\
                          "N: 3";
    ASSERT_STREQ(out.c_str(), expectedString);
    auto loaded = YAML::Load(out.c_str());
    auto emitterStats = loaded.as<Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,false>>();

    ASSERT_TRUE(emitterStats.mean().isApprox(stats.mean()));
    ASSERT_TRUE(emitterStats.standardDeviation().isApprox(stats.standardDeviation()));
    ASSERT_TRUE(emitterStats.cwiseMin().isApprox(stats.cwiseMin()));
    ASSERT_TRUE(emitterStats.cwiseMax().isApprox(stats.cwiseMax()));
    ASSERT_EQ(emitterStats.getTotalWeight(), stats.getTotalWeight());


    YAML::Node node;
    node = stats;
    auto nodeStats = node.as<Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,false>>();

    ASSERT_TRUE(nodeStats.mean().isApprox(stats.mean()));
    ASSERT_TRUE(nodeStats.standardDeviation().isApprox(stats.standardDeviation()));
    ASSERT_TRUE(nodeStats.cwiseMin().isApprox(stats.cwiseMin()));
    ASSERT_TRUE(nodeStats.cwiseMax().isApprox(stats.cwiseMax()));
    ASSERT_EQ(nodeStats.getTotalWeight(), stats.getTotalWeight());
}

TEST_F(AStatisticsTest, YAMLVectorConversion){
    Statistics::RunningStatistics<Eigen::VectorXd,unsigned,false> stats;
    stats.add(vec1);
    stats.add(vec2);
    stats.add(vec3);

    YAML::Emitter out;
    out << stats;

    auto expectedString = "0: [4, 1.154700538379252, 2, 6]\n"\
                          "1: [8, 2.309401076758503, 4, 12]\n"\
                          "2: [12, 3.464101615137755, 6, 18]\n"\
                          "3: [16, 4.618802153517007, 8, 24]\n"\
                          "N: 3";
    ASSERT_STREQ(out.c_str(), expectedString);
    auto loaded = YAML::Load(out.c_str());
    auto emitterStats = loaded.as<Statistics::RunningStatistics<Eigen::VectorXd,unsigned,false>>();

    ASSERT_TRUE(emitterStats.mean().isApprox(stats.mean()));
    ASSERT_TRUE(emitterStats.standardDeviation().isApprox(stats.standardDeviation()));
    ASSERT_TRUE(emitterStats.cwiseMin().isApprox(stats.cwiseMin()));
    ASSERT_TRUE(emitterStats.cwiseMax().isApprox(stats.cwiseMax()));
    ASSERT_EQ(emitterStats.getTotalWeight(), stats.getTotalWeight());


    YAML::Node node;
    node = stats;
    auto nodeStats = node.as<Statistics::RunningStatistics<Eigen::VectorXd,unsigned,false>>();

    ASSERT_TRUE(nodeStats.mean().isApprox(stats.mean()));
    ASSERT_TRUE(nodeStats.standardDeviation().isApprox(stats.standardDeviation()));
    ASSERT_TRUE(nodeStats.cwiseMin().isApprox(stats.cwiseMin()));
    ASSERT_TRUE(nodeStats.cwiseMax().isApprox(stats.cwiseMax()));
    ASSERT_EQ(nodeStats.getTotalWeight(), stats.getTotalWeight());
}

TEST_F(AStatisticsTest, YAMLConversionTriangularExport){
    Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true> stats;
    stats.add(mat1);
    stats.add(mat2);
    stats.add(mat3);

    Eigen::MatrixXd expectedMean(2,2),expectedStandardDeviation(2,2),expectedCwiseMin(2,2),expectedCwiseMax(2,2);
    expectedMean << 0,8,0,0;
    expectedStandardDeviation << 0,4,0,0;
    expectedCwiseMin << 0,4,0,0;
    expectedCwiseMax << 0,12,0,0;

    YAML::Emitter out;
    out << stats;
    auto expectedString = "0:\n"\
                          "  1: [8, 2.309401076758503, 4, 12]\n"\
                          "N: 3";
    ASSERT_STREQ(out.c_str(), expectedString);
    auto loaded = YAML::Load(out.c_str());
    auto emitterStats = loaded.as<Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true>>();

    ASSERT_TRUE(emitterStats.mean().isApprox(expectedMean));
    ASSERT_TRUE(emitterStats.standardDeviation().isApprox(expectedStandardDeviation));
    ASSERT_TRUE(emitterStats.cwiseMin().isApprox(expectedCwiseMin));
    ASSERT_TRUE(emitterStats.cwiseMax().isApprox(expectedCwiseMax));
    ASSERT_EQ(emitterStats.getTotalWeight(), stats.getTotalWeight());


    YAML::Node node;
    node = stats;
    auto nodeStats = node.as<Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true>>();

    ASSERT_TRUE(nodeStats.mean().isApprox(expectedMean));
    ASSERT_TRUE(nodeStats.standardDeviation().isApprox(expectedStandardDeviation));
    ASSERT_TRUE(nodeStats.cwiseMin().isApprox(expectedCwiseMin));
    ASSERT_TRUE(nodeStats.cwiseMax().isApprox(expectedCwiseMax));
    ASSERT_EQ(nodeStats.getTotalWeight(), stats.getTotalWeight());
}
