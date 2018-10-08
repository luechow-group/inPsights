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

    void SetUp() override {
        auto rng = std::default_random_engine {};

        maxima = {
                {1.00, ElectronsVector({{Spin::alpha,{0,0,1.00}},{Spin::beta,{0,0,0}}}),0},
                {1.01, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.01}}}),1},
                {1.02, ElectronsVector({{Spin::alpha,{0,0,1.02}},{Spin::beta,{0,0,0}}}),2},
                {1.03, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.03}}}),3},
                {1.10, ElectronsVector({{Spin::alpha,{0,0,1.10}},{Spin::beta,{0,0,0}}}),4},
                {1.11, ElectronsVector({{Spin::alpha,{0,0,1.11}},{Spin::beta,{0,0,0}}}),5},
                {1.12, ElectronsVector({{Spin::alpha,{0,0,1.12}},{Spin::beta,{0,0,0}}}),6},
                {1.13, ElectronsVector({{Spin::alpha,{0,0,1.13}},{Spin::beta,{0,0,0}}}),7}
        };
    }
};

TEST_F(AGlobalSimilaritySorterTest, OneList) {
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(maxima, similarReferencesVector, 1);
    globalSimilaritySorter.sort();

    ASSERT_EQ((*similarReferencesVector.at(0).repRefIt_).id_,0);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(0).it_).id_,1);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(1).it_).id_,2);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(2).it_).id_,3);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(3).it_).id_,4);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(4).it_).id_,5);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(5).it_).id_,6);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(6).it_).id_,7);
}


TEST_F(AGlobalSimilaritySorterTest, TwoLists) {
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(maxima, similarReferencesVector, 0.1);
    globalSimilaritySorter.sort();

    ASSERT_EQ((*similarReferencesVector.at(0).repRefIt_).id_,0);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(0).it_).id_,1);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(1).it_).id_,2);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(2).it_).id_,3);

    ASSERT_EQ((*similarReferencesVector.at(1).repRefIt_).id_,4);
    ASSERT_EQ((*similarReferencesVector.at(1).similarReferences_.at(0).it_).id_,5);
    ASSERT_EQ((*similarReferencesVector.at(1).similarReferences_.at(1).it_).id_,6);
    ASSERT_EQ((*similarReferencesVector.at(1).similarReferences_.at(2).it_).id_,7);
}

TEST_F(AGlobalSimilaritySorterTest, TwoListsIncrementBorderCase) {
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(maxima, similarReferencesVector, 0.02);
    globalSimilaritySorter.sort();

    ASSERT_EQ((*similarReferencesVector.at(0).repRefIt_).id_,0);
    ASSERT_EQ((*similarReferencesVector.at(0).similarReferences_.at(0).it_).id_,1);

    ASSERT_EQ((*similarReferencesVector.at(1).repRefIt_).id_,2);
    ASSERT_EQ((*similarReferencesVector.at(1).similarReferences_.at(0).it_).id_,3);

    ASSERT_EQ((*similarReferencesVector.at(2).repRefIt_).id_,4);
    ASSERT_EQ((*similarReferencesVector.at(2).similarReferences_.at(0).it_).id_,5);

    ASSERT_EQ((*similarReferencesVector.at(3).repRefIt_).id_,6);
    ASSERT_EQ((*similarReferencesVector.at(3).similarReferences_.at(0).it_).id_,7);
}
