//
// Created by heuer on 18.04.19.
//

#include <gmock/gmock.h>
#include <BestMatchSOAPSimilarityClusterer.h>
#include <BestMatchDistanceSimilarityClusterer.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <Group.h>
#include <Reference.h>
#include <BestMatchSimilarity.h>

using namespace testing;
using namespace SOAP;

class ABestMatchSOAPSimilarityClustererTest : public ::testing::Test {
public:
    Group A, B, C, D, E, F;
    std::vector<Sample> samples;
    Eigen::VectorXd ekin;

    Electron ea,eb;
    AtomsVector atoms;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);

        double a = 1.0, b = a/2.0, c = a*2.0;
        atoms = AtomsVector({{Elements::ElementType::H,{0, 0,-a}},
                             {Elements::ElementType::H,{0, 0,+a}}});

        A = Group({1.0, ElectronsVector({{Spin::alpha, {0, 0,-a}},
                                         {Spin::alpha, {0, 0, a}},
                                         {Spin::beta,  {0, b, 0}}}), 0}); // covalent
        B = Group({1.1, ElectronsVector({{Spin::beta,  {0, b, 0}},
                                         {Spin::alpha, {0, 0, a}},
                                         {Spin::alpha, {0, 0,-a}}}), 1}); // A permuted
        C = Group({1.2, ElectronsVector({{Spin::alpha, {0, 0,-a}},
                                         {Spin::alpha, {0, 0, a}},
                                         {Spin::beta,  {0,-b, 0}}}), 2}); // A reflected
        D = Group({1.3, ElectronsVector({{Spin::alpha, {0, 0,-a}},
                                         {Spin::beta,  {0,-b, 0}},
                                         {Spin::alpha, {0, 0, a}}}), 3}); // A reflected and permuted
        E = Group({1.4, ElectronsVector({{Spin::alpha, {0, 0,-a}},
                                         {Spin::beta,  {0, c, 0}},
                                         {Spin::alpha, {0, 0, a}}}), 4}); // ionic
        F = Group({1.5, ElectronsVector({{Spin::alpha, {0, 0,-a}},
                                         {Spin::beta,  {0, c, 0}},
                                         {Spin::alpha, {0, 0, a}}}), 5}); // E reflected
        ekin = Eigen::VectorXd::Zero(3);

        ea = {Spin::alpha,{0,0,0}};
        eb = {Spin::beta, {0,0,0}};

        samples = {
                {ElectronsVector({ea,ea,eb}), ekin},
                {ElectronsVector({ea,ea,eb}), ekin},
                {ElectronsVector({ea,ea,eb}), ekin},
                {ElectronsVector({ea,ea,eb}), ekin},
                {ElectronsVector({ea,ea,eb}), ekin},
                {ElectronsVector({ea,ea,eb}), ekin}};
    }
};

TEST_F(ABestMatchSOAPSimilarityClustererTest, VerifyTestCluster) {
    ParticleKit::create(atoms, A.representative()->maximum());
    BestMatchSOAPSimilarityClusterer bestMatchSOAPSimilarityClusterer(atoms, samples);

    General::settings.mode = General::Mode::chemical;
    BestMatchSOAPSimilarityClusterer::settings.threshold = 1.0;
    BestMatchDistanceSimilarityClusterer::settings.similarityRadius = 0.1;
    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;

    //General::settings.mode = General::Mode::alchemical;
    //General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;

    auto specA = MolecularSpectrum({atoms, A.representative()->maximum()});
    auto specB = MolecularSpectrum({atoms, B.representative()->maximum()});
    auto specC = MolecularSpectrum({atoms, C.representative()->maximum()});
    auto specD = MolecularSpectrum({atoms, D.representative()->maximum()});
    auto specE = MolecularSpectrum({atoms, E.representative()->maximum()});
    auto specF = MolecularSpectrum({atoms, F.representative()->maximum()});

    auto resAB = BestMatch::SOAPSimilarity::compare(specA,specB,
            BestMatchSOAPSimilarityClusterer::settings.threshold.get(),0.1 );

    ASSERT_EQ(BestMatch::SOAPSimilarity::compare(
            specA,specB,BestMatchSOAPSimilarityClusterer::settings.threshold.get(),0.1 ).metric,1);
    ASSERT_EQ(BestMatch::SOAPSimilarity::compare(
            specA,specC,BestMatchSOAPSimilarityClusterer::settings.threshold.get(),0.1 ).metric,1);
    ASSERT_EQ(BestMatch::SOAPSimilarity::compare(
            specA,specD,BestMatchSOAPSimilarityClusterer::settings.threshold.get(),0.1 ).metric,1);

    ASSERT_LT(BestMatch::SOAPSimilarity::compare(
            specA,specE,BestMatchSOAPSimilarityClusterer::settings.threshold.get(),0.1 ).metric,1);
    ASSERT_LT(BestMatch::SOAPSimilarity::compare(
            specA,specF,BestMatchSOAPSimilarityClusterer::settings.threshold.get(),0.1 ).metric,1);

    ASSERT_EQ(BestMatch::SOAPSimilarity::compare(
            specE,specF,BestMatchSOAPSimilarityClusterer::settings.threshold.get(),0.1 ).metric,1);
    ASSERT_EQ(BestMatch::SOAPSimilarity::compare(
            specF,specE,BestMatchSOAPSimilarityClusterer::settings.threshold.get(),0.1 ).metric,1);
}

TEST_F(ABestMatchSOAPSimilarityClustererTest, TwoClusters) {
    ParticleKit::create(atoms, A.representative()->maximum());
    BestMatchSOAPSimilarityClusterer bestMatchSOAPSimilarityClusterer(atoms, samples);
    
    General::settings.mode = General::Mode::chemical;
    BestMatchSOAPSimilarityClusterer::settings.threshold =  1.0;
    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;

    Group maxima({A,{B,C},D,E,F});

    bestMatchSOAPSimilarityClusterer.cluster(maxima);

    ASSERT_EQ(maxima.size(), 2);
    ASSERT_EQ(maxima[0].size(), 3);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0,1,2,3));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1,2));
    ASSERT_THAT(maxima[0][1][0].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(maxima[0][1][1].allSampleIds(), ElementsAre(2));
    ASSERT_THAT(maxima[0][2].allSampleIds(), ElementsAre(3));

    ASSERT_EQ(maxima[1].size(), 2);
    ASSERT_THAT(maxima[1].allSampleIds(), ElementsAre(4, 5));
    ASSERT_THAT(maxima[1][0].allSampleIds(), ElementsAre(4));
    ASSERT_THAT(maxima[1][1].allSampleIds(), ElementsAre(5));

    // Group expected({{A,{B,C},D},{E,F}});
}

TEST_F(ABestMatchSOAPSimilarityClustererTest, TwoClusters_Alchemical) {
    ParticleKit::create(atoms, A.representative()->maximum());
    BestMatchSOAPSimilarityClusterer bestMatchSOAPSimilarityClusterer(atoms, samples);

    General::settings.mode = General::Mode::alchemical;
    General::settings.pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;
    BestMatchSOAPSimilarityClusterer::settings.threshold =  1.0;
    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;

    Group maxima({A,{B,C},D,E,F});

    bestMatchSOAPSimilarityClusterer.cluster(maxima);

    ASSERT_EQ(maxima.size(), 2);
    ASSERT_EQ(maxima[0].size(), 3);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0,1,2,3));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1,2));
    ASSERT_THAT(maxima[0][1][0].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(maxima[0][1][1].allSampleIds(), ElementsAre(2));
    ASSERT_THAT(maxima[0][2].allSampleIds(), ElementsAre(3));

    ASSERT_EQ(maxima[1].size(), 2);
    ASSERT_THAT(maxima[1].allSampleIds(), ElementsAre(4, 5));
    ASSERT_THAT(maxima[1][0].allSampleIds(), ElementsAre(4));
    ASSERT_THAT(maxima[1][1].allSampleIds(), ElementsAre(5));

    // Group expected({{A,{B,C},D},{E,F}});
}