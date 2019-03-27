//
// Created by Michael Heuer on 10.09.18.
//

#include <gtest/gtest.h>
#include <DensityBasedScan.h>
#include <sstream>
#include <HungarianHelper.h>

using namespace testing;
using namespace Eigen;

template <typename Scalar, typename VectorType>
Scalar euclideanDistance(const VectorType &p1, const VectorType &p2) {
    return (p1 - p2).norm();
}

class ADensityBasedScanTest : public Test {
public:

    std::vector<Eigen::VectorXf> dataFloat;
    std::vector<Eigen::VectorXd> dataDouble;

    void SetUp() override {
        Eigen::VectorXf vf(1);
        vf << 0.0f;
        Eigen::VectorXd vd(1);
        vd << 0.0;

        for (int k = 0; k < 5; ++k) {
            for (int j = 0; j < 5; ++j) {
                dataFloat.push_back(vf);
                dataDouble.push_back(vd);
                vf[0] += 0.1f;
                vd[0] += 0.1;
            }
            vf[0] += 1;
            vd[0] += 1;
        }
    }
};

TEST_F(ADensityBasedScanTest, Float) {
    DensityBasedScan<float, Eigen::VectorXf, euclideanDistance<float, Eigen::VectorXf>>  dbscan(dataFloat);

    auto nClusters = dbscan.findClusters(0.20001, 5);// TODO WHY?
    ASSERT_EQ(nClusters,5);

    auto result = dbscan.getLabels();

    std::vector<int> expected{
        0,0,0,0,0,
        1,1,1,1,1,
        2,2,2,2,2,
        3,3,3,3,3,
        4,4,4,4,4};

    ASSERT_EQ(result, expected);
}

TEST_F(ADensityBasedScanTest, Double) {
    DensityBasedScan<double, Eigen::VectorXd, euclideanDistance<double, Eigen::VectorXd>>  dbscan(dataDouble);

    auto nClusters = dbscan.findClusters(0.20001, 5);// TODO WHY?
    ASSERT_EQ(nClusters,5);

    auto result = dbscan.getLabels();

    std::vector<int> expected{
            0,0,0,0,0,
            1,1,1,1,1,
            2,2,2,2,2,
            3,3,3,3,3,
            4,4,4,4,4};

    ASSERT_EQ(result, expected);
}

TEST_F(ADensityBasedScanTest, MinSizeTooLarge) {
    DensityBasedScan<float, Eigen::VectorXf, euclideanDistance<float, Eigen::VectorXf>> dbscan(dataFloat);

    auto nClusters = dbscan.findClusters(0.20001, 6);// TODO WHY?
    ASSERT_EQ(nClusters,0);

    auto result = dbscan.getLabels();

    std::vector<int> expected{
            -1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1};

    ASSERT_EQ(result, expected);
}

TEST_F(ADensityBasedScanTest, PredictEps) {
    DensityBasedScan<float, Eigen::VectorXf, euclideanDistance<float, Eigen::VectorXf>>  dbscan(dataFloat);

    auto result = dbscan.predictEps(4); // careful => cluster indices start with 0

    auto min = *std::min_element(result.begin(), result.end());
    auto max = *std::max_element(result.begin(), result.end());

    ASSERT_GE(min,0.19999f);// TODO WHY?
    ASSERT_LE(max,1.1);
}

TEST_F(ADensityBasedScanTest, BestMatchNormDistanceFunction) {
    std::vector<ElectronsVector> data(
            {ElectronsVector({{ Spin::alpha, {0, 1, 2}},
                              { Spin::alpha, {0, 0, 0}}}),
             ElectronsVector({{ Spin::alpha, {0, 0, 0}},
                              { Spin::alpha, {0, 4, 6}}})
            });

    DensityBasedScan<double, ElectronsVector, Metrics::bestMatchNorm<Eigen::Infinity,2>>  dbscan(data);

    std::vector<int32_t> expected{0,0};

    auto nClusters = dbscan.findClusters(5.00001, 1); // TODO WHY?

    ASSERT_EQ(nClusters,1);
    ASSERT_EQ(dbscan.getLabels(),expected);
}