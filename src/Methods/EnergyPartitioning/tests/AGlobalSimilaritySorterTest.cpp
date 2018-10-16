//
// Created by Michael Heuer on 25.09.18.
//

#include <gmock/gmock.h>
#include <Reference.h>
#include <Sample.h>
#include <GlobalSimilaritySorter.h>
#include <algorithm>
#include <random>

using namespace testing;

class AGlobalSimilaritySorterTest : public ::testing::Test {
public:
    std::vector<Reference> maxima;
    std::vector<Sample> samples;
    Eigen::VectorXd ekin;
    void SetUp() override {
        ekin.resize(2);
        ekin[0] = 0;
        ekin[1] = 0;

        maxima = {
                {1.00, ElectronsVector({{Spin::alpha,{0,0,1.00}},{Spin::beta,{0,0,0}}}), 0},
                {1.01, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.01}}}), 1},
                {1.02, ElectronsVector({{Spin::alpha,{0,0,1.02}},{Spin::beta,{0,0,0}}}), 2},
                {1.03, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.03}}}), 3},
                {1.10, ElectronsVector({{Spin::alpha,{0,0,1.10}},{Spin::beta,{0,0,0}}}), 4},
                {1.11, ElectronsVector({{Spin::alpha,{0,0,1.11}},{Spin::beta,{0,0,0}}}), 5},
                {1.12, ElectronsVector({{Spin::alpha,{0,0,1.12}},{Spin::beta,{0,0,0}}}), 6},
                {1.13, ElectronsVector({{Spin::alpha,{0,0,1.13}},{Spin::beta,{0,0,0}}}), 7}
        };
        samples = {
                {ElectronsVector({{Spin::alpha,{0,0,1.00}},{Spin::beta,{0,0,0}}}), ekin},
                {ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.01}}}), ekin},
                {ElectronsVector({{Spin::alpha,{0,0,1.02}},{Spin::beta,{0,0,0}}}), ekin},
                {ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.03}}}), ekin},
                {ElectronsVector({{Spin::alpha,{0,0,1.10}},{Spin::beta,{0,0,0}}}), ekin},
                {ElectronsVector({{Spin::alpha,{0,0,1.11}},{Spin::beta,{0,0,0}}}), ekin},
                {ElectronsVector({{Spin::alpha,{0,0,1.12}},{Spin::beta,{0,0,0}}}), ekin},
                {ElectronsVector({{Spin::alpha,{0,0,1.13}},{Spin::beta,{0,0,0}}}), ekin}
        };
    }
};

TEST_F(AGlobalSimilaritySorterTest, OneList) {
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(samples, maxima, similarReferencesVector, 1);
    globalSimilaritySorter.sort();

    ASSERT_EQ((*similarReferencesVector.at(0).repRefIt_).id_,0);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(0)).id_,1);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(1)).id_,2);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(2)).id_,3);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(3)).id_,4);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(4)).id_,5);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(5)).id_,6);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(6)).id_,7);
}


TEST_F(AGlobalSimilaritySorterTest, TwoLists) {
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(samples, maxima, similarReferencesVector, 0.1);
    globalSimilaritySorter.sort();

    ASSERT_EQ((*similarReferencesVector.at(0).repRefIt_).id_,0);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(0)).id_,1);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(1)).id_,2);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(2)).id_,3);

    ASSERT_EQ((*similarReferencesVector.at(1).repRefIt_).id_,4);
    ASSERT_EQ((*similarReferencesVector.at(1).similarReferences_.at(0)).id_,5);
    ASSERT_EQ((*similarReferencesVector.at(1).similarReferences_.at(1)).id_,6);
    ASSERT_EQ((*similarReferencesVector.at(1).similarReferences_.at(2)).id_,7);
}

TEST_F(AGlobalSimilaritySorterTest, TwoListsIncrementBorderCase) {
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(samples, maxima, similarReferencesVector, 0.02);
    globalSimilaritySorter.sort();

    ASSERT_EQ((*similarReferencesVector.at(0).repRefIt_).id_,0);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(0)).id_,1);

    ASSERT_EQ((*similarReferencesVector.at(1).repRefIt_).id_,2);
    ASSERT_EQ((*similarReferencesVector.at(1).similarReferences_.at(0)).id_,3);

    ASSERT_EQ((*similarReferencesVector.at(2).repRefIt_).id_,4);
    ASSERT_EQ((*similarReferencesVector.at(2).similarReferences_.at(0)).id_,5);

    ASSERT_EQ((*similarReferencesVector.at(3).repRefIt_).id_,6);
    ASSERT_EQ((*similarReferencesVector.at(3).similarReferences_.at(0)).id_,7);
}
