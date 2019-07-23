//
// Created by Michael Heuer on 25.09.18.
//

#include <gmock/gmock.h>
#include <Reference.h>
#include <Sample.h>
#include <DistanceClusterer.h>
#include <IdentityClusterer.h>
#include <algorithm>
#include <random>

using namespace testing;

class ABestMatchDistanceSimilarityClustererTest : public ::testing::Test {
public:
    Group maxima;
    std::vector<Sample> samples;
    Eigen::VectorXd ekin;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        IdentityClusterer::settings.identityRadius = 1e-4; // prevent assert

        ekin.resize(2);
        ekin[0] = 0;
        ekin[1] = 0;

        maxima = {
               Group({1.00, ElectronsVector({{Spin::alpha,{0,0,1.00}},{Spin::beta,{0,0,0}}}), 0}),
               Group({1.01, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.01}}}), 1}),
               Group({1.02, ElectronsVector({{Spin::alpha,{0,0,1.02}},{Spin::beta,{0,0,0}}}), 2}),
               Group({1.03, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.03}}}), 3}),
               Group({1.10, ElectronsVector({{Spin::alpha,{0,0,1.10}},{Spin::beta,{0,0,0}}}), 4}),
               Group({1.11, ElectronsVector({{Spin::alpha,{0,0,1.11}},{Spin::beta,{0,0,0}}}), 5}),
               Group({1.12, ElectronsVector({{Spin::alpha,{0,0,1.12}},{Spin::beta,{0,0,0}}}), 6}),
               Group({1.13, ElectronsVector({{Spin::alpha,{0,0,1.13}},{Spin::beta,{0,0,0}}}), 7})
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

TEST_F(ABestMatchDistanceSimilarityClustererTest, OneList) {
    BestMatchDistanceSimilarityClusterer globalSimilaritySorter(samples);
    BestMatchDistanceSimilarityClusterer::settings.similarityRadius = 1;
    BestMatchDistanceSimilarityClusterer::settings.similarityValueIncrement = 1;
    globalSimilaritySorter.cluster(maxima);

    ASSERT_EQ(maxima.size(), 1);
    ASSERT_THAT(maxima.allSampleIds(), ElementsAre(0,1,2,3,4,5,6,7));

    ASSERT_EQ(maxima[0].size(), 8);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0,1,2,3,4,5,6,7));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(maxima[0][2].allSampleIds(), ElementsAre(2));
    ASSERT_THAT(maxima[0][3].allSampleIds(), ElementsAre(3));
    ASSERT_THAT(maxima[0][4].allSampleIds(), ElementsAre(4));
    ASSERT_THAT(maxima[0][5].allSampleIds(), ElementsAre(5));
    ASSERT_THAT(maxima[0][6].allSampleIds(), ElementsAre(6));
    ASSERT_THAT(maxima[0][7].allSampleIds(), ElementsAre(7));

    ASSERT_EQ(maxima.representative()->ownId(), 0);
    ASSERT_EQ(maxima[0].representative()->ownId(), 0);
    ASSERT_EQ(maxima[0][0].representative()->ownId(), 0);
    ASSERT_EQ(maxima[0][1].representative()->ownId(), 1);
    ASSERT_EQ(maxima[0][2].representative()->ownId(), 2);
    ASSERT_EQ(maxima[0][3].representative()->ownId(), 3);
    ASSERT_EQ(maxima[0][4].representative()->ownId(), 4);
    ASSERT_EQ(maxima[0][5].representative()->ownId(), 5);
    ASSERT_EQ(maxima[0][6].representative()->ownId(), 6);
    ASSERT_EQ(maxima[0][7].representative()->ownId(), 7);
}


TEST_F(ABestMatchDistanceSimilarityClustererTest, TwoLists) {
    BestMatchDistanceSimilarityClusterer globalSimilaritySorter(samples);
    BestMatchDistanceSimilarityClusterer::settings.similarityRadius = 0.1;
    BestMatchDistanceSimilarityClusterer::settings.similarityValueIncrement = 1;
    globalSimilaritySorter.cluster(maxima);

    ASSERT_EQ(maxima.size(), 2);
    ASSERT_THAT(maxima.allSampleIds(), ElementsAre(0,1,2,3,4,5,6,7));

    ASSERT_EQ(maxima[0].size(), 4);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0,1,2,3));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(maxima[0][2].allSampleIds(), ElementsAre(2));
    ASSERT_THAT(maxima[0][3].allSampleIds(), ElementsAre(3));

    ASSERT_EQ(maxima[1].size(), 4);
    ASSERT_THAT(maxima[1].allSampleIds(), ElementsAre(4,5,6,7));
    ASSERT_THAT(maxima[1][0].allSampleIds(), ElementsAre(4));
    ASSERT_THAT(maxima[1][1].allSampleIds(), ElementsAre(5));
    ASSERT_THAT(maxima[1][2].allSampleIds(), ElementsAre(6));
    ASSERT_THAT(maxima[1][3].allSampleIds(), ElementsAre(7));

    ASSERT_EQ(maxima.representative()->ownId(), 0);
    ASSERT_EQ(maxima[0].representative()->ownId(), 0);
    ASSERT_EQ(maxima[0][0].representative()->ownId(), 0);
    ASSERT_EQ(maxima[0][1].representative()->ownId(), 1);
    ASSERT_EQ(maxima[0][2].representative()->ownId(), 2);
    ASSERT_EQ(maxima[0][3].representative()->ownId(), 3);

    ASSERT_EQ(maxima[1].representative()->ownId(), 4);
    ASSERT_EQ(maxima[1][0].representative()->ownId(), 4);
    ASSERT_EQ(maxima[1][1].representative()->ownId(), 5);
    ASSERT_EQ(maxima[1][2].representative()->ownId(), 6);
    ASSERT_EQ(maxima[1][3].representative()->ownId(), 7);
}

TEST_F(ABestMatchDistanceSimilarityClustererTest, DISABLED_TwoListsIncrementBorderCase) {
    BestMatchDistanceSimilarityClusterer globalSimilaritySorter(samples);
    BestMatchDistanceSimilarityClusterer::settings.similarityRadius = 0.02;
    BestMatchDistanceSimilarityClusterer::settings.similarityValueIncrement = 1;
    globalSimilaritySorter.cluster(maxima);

    ASSERT_EQ(maxima.size(), 4);
    ASSERT_THAT(maxima.allSampleIds(), ElementsAre(0,1,2,3,4,5,6,7));

    ASSERT_EQ(maxima[0].size(), 2);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0,1));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1));

    ASSERT_EQ(maxima[1].size(), 2);
    ASSERT_THAT(maxima[1].allSampleIds(), ElementsAre(2,3));
    ASSERT_THAT(maxima[1][0].allSampleIds(), ElementsAre(2));
    ASSERT_THAT(maxima[1][1].allSampleIds(), ElementsAre(3));

    ASSERT_EQ(maxima[2].size(), 2);
    ASSERT_THAT(maxima[2].allSampleIds(), ElementsAre(4,5));
    ASSERT_THAT(maxima[2][0].allSampleIds(), ElementsAre(4));
    ASSERT_THAT(maxima[2][1].allSampleIds(), ElementsAre(5));

    ASSERT_EQ(maxima[3].size(), 2);
    ASSERT_THAT(maxima[3].allSampleIds(), ElementsAre(6,7));
    ASSERT_THAT(maxima[3][0].allSampleIds(), ElementsAre(6));
    ASSERT_THAT(maxima[3][1].allSampleIds(), ElementsAre(7));


    ASSERT_EQ(maxima.representative()->ownId(), 0);
    ASSERT_EQ(maxima[0].representative()->ownId(), 0);
    ASSERT_EQ(maxima[0][0].representative()->ownId(), 0);
    ASSERT_EQ(maxima[0][1].representative()->ownId(), 1);

    ASSERT_EQ(maxima[1].representative()->ownId(), 2);
    ASSERT_EQ(maxima[1][0].representative()->ownId(), 2);
    ASSERT_EQ(maxima[1][1].representative()->ownId(), 3);

    ASSERT_EQ(maxima[2].representative()->ownId(), 4);
    ASSERT_EQ(maxima[2][0].representative()->ownId(), 4);
    ASSERT_EQ(maxima[2][1].representative()->ownId(), 5);

    ASSERT_EQ(maxima[3].representative()->ownId(), 6);
    ASSERT_EQ(maxima[3][0].representative()->ownId(), 6);
    ASSERT_EQ(maxima[3][1].representative()->ownId(), 7);

}
