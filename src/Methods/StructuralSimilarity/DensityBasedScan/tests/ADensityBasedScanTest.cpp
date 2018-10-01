//
// Created by Michael Heuer on 10.09.18.
//


#include <gtest/gtest.h>
#include <DensityBasedScan.h>
#include <sstream>

using namespace testing;
using namespace Eigen;

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
    DensityBasedScan<float, Eigen::VectorXf> dbscan(dataFloat);

    auto nClusters = dbscan.findClusters(0.20001, 5);
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
    DensityBasedScan<double, Eigen::VectorXd> dbscan(dataDouble);

    auto nClusters = dbscan.findClusters(0.20001, 5);
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
    DensityBasedScan<float, Eigen::VectorXf> dbscan(dataFloat);

    auto nClusters = dbscan.findClusters(0.20001, 6);
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
    DensityBasedScan<float, Eigen::VectorXf> dbscan(dataFloat);

    auto result = dbscan.predictEps(4); // careful => cluster indices start with 0

    auto min = *std::min_element(result.begin(), result.end());
    auto max = *std::max_element(result.begin(), result.end());

    ASSERT_GE(min,0.19999f); // ???? why
    ASSERT_LE(max,1.1);
}