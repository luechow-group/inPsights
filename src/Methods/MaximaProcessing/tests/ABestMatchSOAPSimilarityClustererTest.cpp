//
// Created by heuer on 18.04.19.
//

#include <gmock/gmock.h>
#include <BestMatchSOAPSimilarityClusterer.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <Group.h>
#include <Reference.h>

using namespace testing;
using namespace SOAP;

class ABestMatchSOAPSimilarityClustererTest : public ::testing::Test {
public:
    Group A, B, C, D, E, F;
    std::vector<Sample> samples;
    Eigen::VectorXd ekin;

    Electron a,b;
    AtomsVector atoms;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);

        atoms = AtomsVector({{Elements::ElementType::H,{0, 0,-0.5}},
                             {Elements::ElementType::H,{0, 0,+0.5}}});

        A = Group({1.0, ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                                         {Spin::alpha, {0, 0, 0.5}},
                                         {Spin::beta,  {0, 0, 0.3}}}), 0});
        B = Group({1.1, ElectronsVector({{Spin::beta,  {0, 0, 0.3}},
                                         {Spin::alpha, {0, 0, 0.5}},
                                         {Spin::alpha, {0, 0,-0.5}}}), 1});
        C = Group({1.2, ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                                         {Spin::alpha, {0, 0, 0.5}},
                                         {Spin::beta,  {0, 0,-0.3}}}), 2});
        D = Group({1.3, ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                                         {Spin::beta,  {0, 0,-0.3}},
                                         {Spin::alpha, {0, 0, 0.5}}}), 3});
        E = Group({1.3, ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                                         {Spin::beta,  {0, 0,-0.5}},
                                         {Spin::alpha, {0, 0, 0.5}}}), 4});
        F = Group({1.3, ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                                         {Spin::beta,  {0, 0, 0.5}},
                                         {Spin::alpha, {0, 0, 0.5}}}), 5});
        ekin = Eigen::VectorXd::Zero(3);

        a = {Spin::alpha,{0,0,0}};
        b = {Spin::beta, {0,0,0}};

        samples = {
                {ElectronsVector({a,a,b}), ekin},
                {ElectronsVector({a,a,b}), ekin},
                {ElectronsVector({a,a,b}), ekin},
                {ElectronsVector({a,a,b}), ekin},
                {ElectronsVector({a,a,b}), ekin},
                {ElectronsVector({a,a,b}), ekin}};
    }
};

TEST_F(ABestMatchSOAPSimilarityClustererTest, TwoClusters) {
    ParticleKit::create(atoms, A.representative()->maximum());
    BestMatchSOAPSimilarityClusterer bestMatchSOAPSimilarityClusterer(atoms, samples);
    
    General::settings.mode = General::Mode::chemical;

    Group maxima({A,{B,C},D,E,F});

    bestMatchSOAPSimilarityClusterer.cluster(maxima);

    ASSERT_EQ(maxima.size(), 2);
    ASSERT_EQ(maxima[0].size(), 4);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0,1,2,3));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(maxima[0][2].allSampleIds(), ElementsAre(2));
    ASSERT_THAT(maxima[0][3].allSampleIds(), ElementsAre(3));

    ASSERT_EQ(maxima[1].size(), 2);
    ASSERT_THAT(maxima[1].allSampleIds(), ElementsAre(4, 5));
    ASSERT_THAT(maxima[1][0].allSampleIds(), ElementsAre(4));
    ASSERT_THAT(maxima[1][1].allSampleIds(), ElementsAre(5));

}

