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
    Clustering::DensityBasedScan<float> dbscan(dataFloat);

    auto nClusters = dbscan.predict(0.20001,5);
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
    Clustering::DensityBasedScan<double> dbscan(dataDouble);

    auto nClusters = dbscan.predict(0.20001,5);
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
    Clustering::DensityBasedScan<float> dbscan(dataFloat);

    auto nClusters = dbscan.predict(0.20001,6);
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
