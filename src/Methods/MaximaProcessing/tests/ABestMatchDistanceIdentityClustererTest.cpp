//
// Created by Michael Heuer on 25.09.18.
//

#include <gmock/gmock.h>
#include <Reference.h>
#include <Sample.h>
#include <BestMatchDistanceIdentityClusterer.h>
#include <BestMatchDistanceSimilarityClusterer.h>
#include <algorithm>
#include <random>

using namespace testing;

class ABestMatchDistanceIdentityClustererTest : public ::testing::Test {
public:
    Group tripletMaxima, singletMaxima;
    std::vector<Sample> tripletSamples, singletSamples;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        BestMatchDistanceSimilarityClusterer::settings.similarityRadius = 10; // prevent assert

        auto rng = std::default_random_engine(static_cast<unsigned long>(std::clock()));

        // Multiplicity = 3: Spin flip is not possible.
        tripletMaxima = {
                Group({1.00, ElectronsVector({{Spin::alpha,{0,0,1.00}},{Spin::alpha,{0,0,0}}}),0}),
                Group({1.01, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,1.01}}}),1}),
                Group({1.02, ElectronsVector({{Spin::alpha,{0,0,1.02}},{Spin::alpha,{0,0,0}}}),2}),
                Group({1.03, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,1.03}}}),3}),
                Group({1.10, ElectronsVector({{Spin::alpha,{0,0,1.10}},{Spin::alpha,{0,0,0}}}),4}),
                Group({1.11, ElectronsVector({{Spin::alpha,{0,0,1.11}},{Spin::alpha,{0,0,0}}}),5}),
                Group({1.12, ElectronsVector({{Spin::alpha,{0,0,1.12}},{Spin::alpha,{0,0,0}}}),6}),
                Group({1.13, ElectronsVector({{Spin::alpha,{0,0,1.13}},{Spin::alpha,{0,0,0}}}),7}),
        };
        std::shuffle(std::begin(tripletMaxima), std::end(tripletMaxima), rng);

        // Multiplicity = 1: Spin flip is possible.
        singletMaxima = {
                Group({1.00, ElectronsVector({{Spin::alpha,{0,0,1.00}},{Spin::beta,{0,0,0}}}),0}),
                Group({1.01, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.01}}}),1}),
                Group({1.02, ElectronsVector({{Spin::alpha,{0,0,1.02}},{Spin::beta,{0,0,0}}}),2}),
                Group({1.03, ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::beta,{0,0,1.03}}}),3}),
                Group({1.10, ElectronsVector({{Spin::alpha,{0,0,1.10}},{Spin::beta,{0,0,0}}}),4}),
                Group({1.11, ElectronsVector({{Spin::alpha,{0,0,1.11}},{Spin::beta,{0,0,0}}}),5}),
                Group({1.12, ElectronsVector({{Spin::alpha,{0,0,1.12}},{Spin::beta,{0,0,0}}}),6}),
                Group({1.13, ElectronsVector({{Spin::alpha,{0,0,1.13}},{Spin::beta,{0,0,0}}}),7}),
        };
        std::shuffle(std::begin(singletMaxima), std::end(singletMaxima), rng);

        for (auto& i : tripletMaxima){
            Sample s(ElectronsVector({{i.representative()->maximum().typesVector()[0],{0,0,0}},
                                      {i.representative()->maximum().typesVector()[1],{0,0,0}}}), Eigen::VectorXd::Random(2));
            tripletSamples.emplace_back(std::move(s));
        }

        for (auto& i : singletMaxima){
            Sample s(ElectronsVector({{i.representative()->maximum().typesVector()[0],{0,0,0}},
                                      {i.representative()->maximum().typesVector()[1],{0,0,0}}}), Eigen::VectorXd::Random(2));
            singletSamples.emplace_back(std::move(s));
        }
    }
};

TEST_F(ABestMatchDistanceIdentityClustererTest, OneListTriplet) {
    BestMatchDistanceIdentityClusterer globalIdentiySorter(tripletSamples);
    BestMatchDistanceIdentityClusterer::settings.identityRadius = 2;
    BestMatchDistanceIdentityClusterer::settings.identityValueIncrement = 1;
    globalIdentiySorter.cluster(tripletMaxima);

    ASSERT_THAT(tripletMaxima.representative()->sampleIds(), ElementsAre(0,1,2,3,4,5,6,7));
}

TEST_F(ABestMatchDistanceIdentityClustererTest, OneListSinglet) {
    BestMatchDistanceIdentityClusterer globalIdentiySorter(singletSamples);
    BestMatchDistanceIdentityClusterer::settings.identityRadius = 2;
    BestMatchDistanceIdentityClusterer::settings.identityValueIncrement = 1;
    globalIdentiySorter.cluster(singletMaxima);

    ASSERT_THAT(singletMaxima.at(0).representative()->sampleIds(), ElementsAre(0,1,2,3,4,5,6,7));
}

TEST_F(ABestMatchDistanceIdentityClustererTest, TwoListsTriplet) {
    BestMatchDistanceIdentityClusterer globalIdentiySorter(tripletSamples);
    BestMatchDistanceIdentityClusterer::settings.identityRadius = 1;
    BestMatchDistanceIdentityClusterer::settings.identityValueIncrement = 0.05;
    globalIdentiySorter.cluster(tripletMaxima);

    ASSERT_THAT(tripletMaxima.at(0).representative()->sampleIds(), ElementsAre(0,1,2,3));
    ASSERT_THAT(tripletMaxima.at(1).representative()->sampleIds(), ElementsAre(4,5,6,7));
}

TEST_F(ABestMatchDistanceIdentityClustererTest, TwoListsSinglet) {
    BestMatchDistanceIdentityClusterer globalIdentiySorter(singletSamples);
    BestMatchDistanceIdentityClusterer::settings.identityRadius = 1;
    BestMatchDistanceIdentityClusterer::settings.identityValueIncrement = 0.05;
    globalIdentiySorter.cluster(singletMaxima);

    ASSERT_THAT(singletMaxima.at(0).representative()->sampleIds(), ElementsAre(0,1,2,3));
    ASSERT_THAT(singletMaxima.at(1).representative()->sampleIds(), ElementsAre(4,5,6,7));
}

TEST_F(ABestMatchDistanceIdentityClustererTest, TwoListsIncrementBorderCaseTriplet) {
    BestMatchDistanceIdentityClusterer globalIdentiySorter(tripletSamples);
    BestMatchDistanceIdentityClusterer::settings.identityRadius = 1;
    BestMatchDistanceIdentityClusterer::settings.identityValueIncrement = 0.1;
    globalIdentiySorter.cluster(tripletMaxima);

    ASSERT_THAT(tripletMaxima.at(0).representative()->sampleIds(), ElementsAre(0,1,2,3,4));
    ASSERT_THAT(tripletMaxima.at(1).representative()->sampleIds(), ElementsAre(5,6,7));
}

TEST_F(ABestMatchDistanceIdentityClustererTest, TwoListsIncrementBorderCaseSinglet) {
    BestMatchDistanceIdentityClusterer globalIdentiySorter(singletSamples);
    BestMatchDistanceIdentityClusterer::settings.identityRadius = 1;
    BestMatchDistanceIdentityClusterer::settings.identityValueIncrement = 0.1;
    globalIdentiySorter.cluster(singletMaxima);

    ASSERT_THAT(singletMaxima.at(0).representative()->sampleIds(), ElementsAre(0,1,2,3,4));
    ASSERT_THAT(singletMaxima.at(1).representative()->sampleIds(), ElementsAre(5,6,7));
}
