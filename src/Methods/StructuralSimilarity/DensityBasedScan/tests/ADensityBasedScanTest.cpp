//
// Created by Michael Heuer on 10.09.18.
//

#include <gtest/gtest.h>
#include <DensityBasedScan.h>
#include <BestMatchDistance.h>
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

    template <typename Scalar, typename VectorType>
    static Scalar euclideanDistance(const VectorType &p1, const VectorType &p2) {
        return (p1 - p2).norm();
    }

    static double bestMatchDistance(const ElectronsVector &e1, const ElectronsVector &e2) {
        return BestMatch::Distance::compare<Eigen::Infinity, 2>(e1, e2).metric;
    };
};

TEST_F(ADensityBasedScanTest, Float) {
    DensityBasedScan<float, Eigen::VectorXf, ADensityBasedScanTest::euclideanDistance<float, Eigen::VectorXf>>  dbscan(dataFloat);

    auto result = dbscan.findClusters(0.20001, 5);// TODO WHY?
    ASSERT_EQ(result.numberOfClusters,5);

    std::vector<int> expected{
        0,0,0,0,0,
        1,1,1,1,1,
        2,2,2,2,2,
        3,3,3,3,3,
        4,4,4,4,4};

    ASSERT_EQ(result.labels, expected);
}

TEST_F(ADensityBasedScanTest, Double) {
    DensityBasedScan<double, Eigen::VectorXd, ADensityBasedScanTest::euclideanDistance<double, Eigen::VectorXd>>  dbscan(dataDouble);

    auto result = dbscan.findClusters(0.20001, 5);// TODO WHY?
    ASSERT_EQ(result.numberOfClusters,5);

    std::vector<int> expected{
            0,0,0,0,0,
            1,1,1,1,1,
            2,2,2,2,2,
            3,3,3,3,3,
            4,4,4,4,4};


    ASSERT_EQ(result.labels, expected);
}

TEST_F(ADensityBasedScanTest, MinSizeTooLarge) {
    DensityBasedScan<float, Eigen::VectorXf, ADensityBasedScanTest::euclideanDistance<float, Eigen::VectorXf>> dbscan(dataFloat);

    auto result = dbscan.findClusters(0.20001, 6);// TODO WHY?
    ASSERT_EQ(result.numberOfClusters, 0);
    
    std::vector<int> expected{
            -1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1,
            -1,-1,-1,-1,-1};

    ASSERT_EQ(result.labels, expected);
}

TEST_F(ADensityBasedScanTest, PredictEps) {
    DensityBasedScan<float, Eigen::VectorXf, ADensityBasedScanTest::euclideanDistance<float, Eigen::VectorXf>>  dbscan(dataFloat);

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

    DensityBasedScan<double, ElectronsVector, ADensityBasedScanTest::bestMatchDistance>  dbscan(data);

    std::vector<int32_t> expected{0,0};

    auto result = dbscan.findClusters(5.00001, 1); // TODO WHY?

    ASSERT_EQ(result.numberOfClusters,1);
    ASSERT_EQ(result.labels,expected);
}