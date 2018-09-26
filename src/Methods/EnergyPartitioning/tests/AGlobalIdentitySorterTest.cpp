//
// Created by Michael Heuer on 25.09.18.
//

#include <gmock/gmock.h>
#include <Reference.h>
#include <Sample.h>
#include <GlobalIdentitySorter.h>

using namespace testing;

class AGlobalIdentitySorterTest : public ::testing::Test {
public:
    std::vector<Reference> system1;
    std::vector<Sample> samples;

    void SetUp() override {

        // multiplicity = 2, no spin flip possible
        system1 = { 
                {1.00,ElectronsVector({{Spin::alpha,{0,0,1.00}}}),0},
                {1.01,ElectronsVector({{Spin::alpha,{0,0,1.01}}}),1},
                {1.02,ElectronsVector({{Spin::alpha,{0,0,1.02}}}),2},
                {1.03,ElectronsVector({{Spin::alpha,{0,0,1.03}}}),3},
                {1.10,ElectronsVector({{Spin::alpha,{0,0,1.10}}}),4},
                {1.11,ElectronsVector({{Spin::alpha,{0,0,1.11}}}),5},
                {1.12,ElectronsVector({{Spin::alpha,{0,0,1.12}}}),6},
                {1.13,ElectronsVector({{Spin::alpha,{0,0,1.13}}}),7},
        };

        for (auto& i : system1){
            Sample s(ElectronsVector({{Spin::alpha,{0,0,0}}}), Eigen::VectorXd::Random(1));
            samples.emplace_back(std::move(s));
        }
    }
};

TEST_F(AGlobalIdentitySorterTest, OneList) {
    GlobalIdentiySorter globalIdentiySorter(system1, samples, 1, 2);
    globalIdentiySorter.sort();

    ASSERT_THAT(system1.at(0).associatedSampleIds_, ElementsAre(1,2,3,4,5,6,7));
}

TEST_F(AGlobalIdentitySorterTest, TwoLists) {
    GlobalIdentiySorter globalIdentiySorter(system1, samples, 0.05, 1);
    globalIdentiySorter.sort();

    ASSERT_THAT(system1.at(0).associatedSampleIds_, ElementsAre(1,2,3));
    ASSERT_THAT(system1.at(1).associatedSampleIds_, ElementsAre(5,6,7));
}

TEST_F(AGlobalIdentitySorterTest, TwoListsIncrementBorderCase) {
    GlobalIdentiySorter globalIdentiySorter(system1, samples, 0.1, 1);
    globalIdentiySorter.sort();

    ASSERT_THAT(system1.at(0).associatedSampleIds_, ElementsAre(1,2,3,4));
    ASSERT_THAT(system1.at(1).associatedSampleIds_, ElementsAre(6,7));
}
