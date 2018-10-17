//
// Created by Michael Heuer on 25.09.18.
//

#include <gmock/gmock.h>
#include <Reference.h>
#include <Sample.h>
#include <GlobalIdentitySorter.h>
#include <algorithm>
#include <random>

using namespace testing;

class AGlobalIdentitySorterTest : public ::testing::Test {
public:
    std::vector<Reference> tripletMaxima, singletMaxima;
    std::vector<Sample> tripletSamples, singletSamples;

    void SetUp() override {
        auto rng = std::default_random_engine {};

        // Multiplicity = 3: Spin flip is not possible.
        tripletMaxima = {
                {1.00, ElectronsVector({{Spin::alpha,{0,0,1.00}},{Spin::alpha,{0,0,0}}}),0},
                {1.01, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,1.01}}}),1},
                {1.02, ElectronsVector({{Spin::alpha,{0,0,1.02}},{Spin::alpha,{0,0,0}}}),2},
                {1.03, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,1.03}}}),3},
                {1.10, ElectronsVector({{Spin::alpha,{0,0,1.10}},{Spin::alpha,{0,0,0}}}),4},
                {1.11, ElectronsVector({{Spin::alpha,{0,0,1.11}},{Spin::alpha,{0,0,0}}}),5},
                {1.12, ElectronsVector({{Spin::alpha,{0,0,1.12}},{Spin::alpha,{0,0,0}}}),6},
                {1.13, ElectronsVector({{Spin::alpha,{0,0,1.13}},{Spin::alpha,{0,0,0}}}),7},
        };
        std::shuffle(std::begin(tripletMaxima), std::end(tripletMaxima), rng);

        // Multiplicity = 1: Spin flip is possible.
        singletMaxima = {
                {1.00, ElectronsVector({{Spin::alpha,{0,0,1.00}},{Spin::beta,{0,0,0}}}),0},
                {1.01, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.01}}}),1},
                {1.02, ElectronsVector({{Spin::alpha,{0,0,1.02}},{Spin::beta,{0,0,0}}}),2},
                {1.03, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.03}}}),3},
                {1.10, ElectronsVector({{Spin::alpha,{0,0,1.10}},{Spin::beta,{0,0,0}}}),4},
                {1.11, ElectronsVector({{Spin::alpha,{0,0,1.11}},{Spin::beta,{0,0,0}}}),5},
                {1.12, ElectronsVector({{Spin::alpha,{0,0,1.12}},{Spin::beta,{0,0,0}}}),6},
                {1.13, ElectronsVector({{Spin::alpha,{0,0,1.13}},{Spin::beta,{0,0,0}}}),7},
        };
        std::shuffle(std::begin(singletMaxima), std::end(singletMaxima), rng);

        for (auto& i : tripletMaxima){
            Sample s(ElectronsVector({{i.maximum_.typesVector()[0],{0,0,0}},
                                      {i.maximum_.typesVector()[1],{0,0,0}}}), Eigen::VectorXd::Random(2));
            tripletSamples.emplace_back(std::move(s));
        }

        for (auto& i : singletMaxima){
            Sample s(ElectronsVector({{i.maximum_.typesVector()[0],{0,0,0}},
                                      {i.maximum_.typesVector()[1],{0,0,0}}}), Eigen::VectorXd::Random(2));
            singletSamples.emplace_back(std::move(s));
        }
    }
};

TEST_F(AGlobalIdentitySorterTest, OneListTriplet) {
    GlobalIdentiySorter globalIdentiySorter(tripletMaxima, tripletSamples, 1, 2);
    globalIdentiySorter.sort();

    ASSERT_THAT(tripletMaxima.at(0).sampleIds_, ElementsAre(0,1,2,3,4,5,6,7));
}

TEST_F(AGlobalIdentitySorterTest, OneListSinglet) {
    GlobalIdentiySorter globalIdentiySorter(singletMaxima, singletSamples, 1, 2);
    globalIdentiySorter.sort();

    ASSERT_THAT(singletMaxima.at(0).sampleIds_, ElementsAre(0,1,2,3,4,5,6,7));
}

TEST_F(AGlobalIdentitySorterTest, TwoListsTriplet) {
    GlobalIdentiySorter globalIdentiySorter(tripletMaxima, tripletSamples, 0.05, 1);
    globalIdentiySorter.sort();

    ASSERT_THAT(tripletMaxima.at(0).sampleIds_, ElementsAre(0,1,2,3));
    ASSERT_THAT(tripletMaxima.at(1).sampleIds_, ElementsAre(4,5,6,7));
}

TEST_F(AGlobalIdentitySorterTest, TwoListsSinglet) {
    GlobalIdentiySorter globalIdentiySorter(singletMaxima, singletSamples, 0.05, 1);
    globalIdentiySorter.sort();

    ASSERT_THAT(singletMaxima.at(0).sampleIds_, ElementsAre(0,1,2,3));
    ASSERT_THAT(singletMaxima.at(1).sampleIds_, ElementsAre(4,5,6,7));
}

TEST_F(AGlobalIdentitySorterTest, TwoListsIncrementBorderCaseTriplet) {
    GlobalIdentiySorter globalIdentiySorter(tripletMaxima, tripletSamples, 0.1, 1);
    globalIdentiySorter.sort();

    ASSERT_THAT(tripletMaxima.at(0).sampleIds_, ElementsAre(0,1,2,3,4));
    ASSERT_THAT(tripletMaxima.at(1).sampleIds_, ElementsAre(5,6,7));
}

TEST_F(AGlobalIdentitySorterTest, TwoListsIncrementBorderCaseSinglet) {
    GlobalIdentiySorter globalIdentiySorter(singletMaxima, singletSamples, 0.1, 1);
    globalIdentiySorter.sort();

    ASSERT_THAT(singletMaxima.at(0).sampleIds_, ElementsAre(0,1,2,3,4));
    ASSERT_THAT(singletMaxima.at(1).sampleIds_, ElementsAre(5,6,7));
}
