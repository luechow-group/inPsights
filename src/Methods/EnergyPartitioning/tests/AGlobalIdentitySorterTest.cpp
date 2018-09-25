//
// Created by Michael Heuer on 25.09.18.
//

#include <gtest/gtest.h>
#include <Reference.h>
#include <Sample.h>
#include <GlobalIdentitySorter.h>

using namespace testing;

class AGlobalIdentitySorterTest : public ::testing::Test {
public:
    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;

    void SetUp() override {

        globallyIdenticalMaxima = {
                {1.00,ElectronsVector({{Spin::alpha,{0,0,1.00}}}),0},
                {1.01,ElectronsVector({{Spin::alpha,{0,0,1.01}}}),1},
                {1.02,ElectronsVector({{Spin::alpha,{0,0,1.02}}}),2},
                {1.03,ElectronsVector({{Spin::alpha,{0,0,1.03}}}),3},
                {1.10,ElectronsVector({{Spin::alpha,{0,0,1.10}}}),4},
                {1.11,ElectronsVector({{Spin::alpha,{0,0,1.11}}}),5},
                {1.12,ElectronsVector({{Spin::alpha,{0,0,1.12}}}),6},
                {1.13,ElectronsVector({{Spin::alpha,{0,0,1.13}}}),7},
        };

        for (auto& i : globallyIdenticalMaxima){
            Sample s(ElectronsVector({{Spin::alpha,{0,0,0}}}), Eigen::VectorXd::Random(1));

            samples.emplace_back(std::move(s));
        }
    }
};

TEST_F(AGlobalIdentitySorterTest, sort) {
    std::cout << globallyIdenticalMaxima.size() << std::endl;

    double increment = 0.1;
    double identicalDistThresh = 1;
    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, increment ,identicalDistThresh);

    globalIdentiySorter.sort();

    std::cout << globallyIdenticalMaxima.size() << std::endl;

}