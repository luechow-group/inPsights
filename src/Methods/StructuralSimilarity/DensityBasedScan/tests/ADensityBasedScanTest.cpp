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
    auto dataPtr = std::make_shared<Clustering::Dataset<float>>(dataFloat);
    Clustering::DensityBasedScan<float> dbscan(dataPtr);
    dbscan.fit();

    dbscan.predict(0.20001,4);
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
    auto dataPtr = std::make_shared<Clustering::Dataset<float>>(dataFloat);
    Clustering::DensityBasedScan<float> dbscan(dataPtr);
    dbscan.fit();

    dbscan.predict(0.20001,4);
    auto result = dbscan.getLabels();

    std::vector<int> expected{
            0,0,0,0,0,
            1,1,1,1,1,
            2,2,2,2,2,
            3,3,3,3,3,
            4,4,4,4,4};

    ASSERT_EQ(result, expected);
}
