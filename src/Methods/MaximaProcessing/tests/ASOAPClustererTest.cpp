/* Copyright (C) 2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <gmock/gmock.h>
#include <SOAPClusterer.h>
#include <PreClusterer.h>
#include <SOAPSettings.h>
#include <TestMolecules.h>
#include <Group.h>
#include <Reference.h>
#include <BestMatchSimilarity.h>

using namespace testing;
using namespace SOAP;

class ASOAPClustererTest : public ::testing::Test {
public:
    Group A, B, C, D, E, F;
    std::vector<Sample> samples;
    Eigen::VectorXd ekin;
    double eps = SOAP::General::settings.comparisonEpsilon();

    Electron ea, eb;
    AtomsVector atoms;

    void SetUp() override {
        //spdlog::set_level(spdlog::level::debug);

        double a = 1.0, b = a / 2.0, c = a * 2.0;
        atoms = AtomsVector({
            {Elements::ElementType::H, {0, 0, -a}},
            {Elements::ElementType::H, {0, 0, +a}}});

        A = Group({1.0, ElectronsVector({
            {Spin::alpha, {0, 0,-a}},
            {Spin::alpha, {0, 0, a}},
            {Spin::beta,  {0, b, 0}}}), 0}); // covalent
        B = Group({1.1, ElectronsVector({
            {Spin::beta,  {0, b, 0}},
            {Spin::alpha, {0, 0, a-0.01}},
            {Spin::alpha, {0, 0,-a+0.01}}}), 1}); // A permuted and shifted
        C = Group({1.1, ElectronsVector({
            {Spin::beta,  {0, b, 0}},
            {Spin::alpha, {0, 0, a+0.01}},
            {Spin::alpha, {0, 0,-a-0.01}}}), 2}); // A permuted and shifted
        D = Group({1.3, ElectronsVector({
            {Spin::alpha, {0, 0,-a}},
            {Spin::beta,  {0,-b, 0}},
            {Spin::alpha, {0, 0, a}}}), 3}); // A reflected and permuted
        E = Group({1.4, ElectronsVector({
            {Spin::alpha, {0, 0,-a}},
            {Spin::beta,  {0, c, 0}},
            {Spin::alpha, {0, 0, a}}}), 4}); // ionic
        F = Group({1.5, ElectronsVector({
            {Spin::alpha, {0, 0,-a}},
            {Spin::beta,  {0, c, 0}},
            {Spin::alpha, {0, 0, a}}}), 5}); // E reflected
        ekin = Eigen::VectorXd::Zero(3);

        ea = {Spin::alpha, {0, 0, 0}};
        eb = {Spin::beta, {0, 0, 0}};

        samples = {
                {ElectronsVector({ea, ea, eb}), ekin},
                {ElectronsVector({ea, ea, eb}), ekin},
                {ElectronsVector({ea, ea, eb}), ekin},
                {ElectronsVector({ea, ea, eb}), ekin},
                {ElectronsVector({ea, ea, eb}), ekin},
                {ElectronsVector({ea, ea, eb}), ekin}};
    }
};

TEST_F(ASOAPClustererTest, VerifyTestCluster) {
    ParticleKit::create(atoms, A.representative()->maximum());

    General::settings.mode = General::Mode::chemical;
    double soapThreshold = 0.99;
    double distanceTolerance = 0.1;

    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;

    auto specA = MolecularSpectrum({atoms, A.representative()->maximum()});
    auto specB = MolecularSpectrum({atoms, B.representative()->maximum()});
    auto specC = MolecularSpectrum({atoms, C.representative()->maximum()});
    auto specD = MolecularSpectrum({atoms, D.representative()->maximum()});
    auto specE = MolecularSpectrum({atoms, E.representative()->maximum()});
    auto specF = MolecularSpectrum({atoms, F.representative()->maximum()});
    
    ASSERT_LT(BestMatch::SOAPSimilarity::compare(
            specA, specB, soapThreshold, distanceTolerance, eps).metric, 1);

    ASSERT_LT(BestMatch::SOAPSimilarity::compare(
            specA, specC, soapThreshold, distanceTolerance, eps).metric, 1);

    ASSERT_EQ(BestMatch::SOAPSimilarity::compare(
            specA, specD, soapThreshold, distanceTolerance, eps).metric, 1);

    ASSERT_LT(BestMatch::SOAPSimilarity::compare(
            specA, specE, soapThreshold, distanceTolerance, eps).metric, 1);
    ASSERT_LT(BestMatch::SOAPSimilarity::compare(
            specA, specF, soapThreshold, distanceTolerance, eps).metric, 1);

    ASSERT_EQ(BestMatch::SOAPSimilarity::compare(
            specE, specF, soapThreshold, distanceTolerance, eps).metric, 1);
    ASSERT_EQ(BestMatch::SOAPSimilarity::compare(
            specF, specE, soapThreshold, distanceTolerance, eps).metric, 1);
}

TEST_F(ASOAPClustererTest, TwoClusters) {
    ParticleKit::create(atoms, A.representative()->maximum());
    SOAPClusterer clusterer(atoms, samples);

    General::settings.mode = General::Mode::chemical;
    SOAPClusterer::settings.similarityThreshold = 0.97;
    SOAPClusterer::settings.distanceMatrixCovarianceTolerance = 0.1;
    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;

    Group maxima({A, {B, C}, D, E, F});

    clusterer.cluster(maxima);

    ASSERT_EQ(maxima.size(), 2);
    ASSERT_EQ(maxima[0].size(), 3);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0, 1, 2, 3));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1, 2));
    ASSERT_THAT(maxima[0][1][0].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(maxima[0][1][1].allSampleIds(), ElementsAre(2));
    ASSERT_THAT(maxima[0][2].allSampleIds(), ElementsAre(3));

    ASSERT_EQ(maxima[1].size(), 2);
    ASSERT_THAT(maxima[1].allSampleIds(), ElementsAre(4, 5));
    ASSERT_THAT(maxima[1][0].allSampleIds(), ElementsAre(4));
    ASSERT_THAT(maxima[1][1].allSampleIds(), ElementsAre(5));

    // Group expected({{A,{B,C},D},{E,F}});
}

TEST_F(ASOAPClustererTest, TwoClusters_Alchemical) {
    ParticleKit::create(atoms, A.representative()->maximum());
    SOAPClusterer clusterer(atoms, samples);

    General::settings.mode = General::Mode::alchemical;
    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    SOAPClusterer::settings.similarityThreshold = 0.97;
    SOAPClusterer::settings.distanceMatrixCovarianceTolerance = 0.1;
    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;

    Group maxima({A, {B, C}, D, E, F});

    clusterer.cluster(maxima);

    ASSERT_EQ(maxima.size(), 2);
    ASSERT_EQ(maxima[0].size(), 3);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0, 1, 2, 3));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1, 2));
    ASSERT_THAT(maxima[0][1][0].allSampleIds(), ElementsAre(1));
    ASSERT_THAT(maxima[0][1][1].allSampleIds(), ElementsAre(2));
    ASSERT_THAT(maxima[0][2].allSampleIds(), ElementsAre(3));

    ASSERT_EQ(maxima[1].size(), 2);
    ASSERT_THAT(maxima[1].allSampleIds(), ElementsAre(4, 5));
    ASSERT_THAT(maxima[1][0].allSampleIds(), ElementsAre(4));
    ASSERT_THAT(maxima[1][1].allSampleIds(), ElementsAre(5));

    // Group expected({{A,{B,C},D},{E,F}});
}


TEST_F(ASOAPClustererTest, BH3_Alchemical) {
    using namespace TestMolecules;
    auto nuclei = BH3::nuclei.atoms();

    const MolecularGeometry B = {
            nuclei,
            ElectronsVector({
                                    {Spin::alpha, nuclei.positionsVector()[1]},
                                    {Spin::beta, nuclei.positionsVector()[2]},
                                    {Spin::beta,inbetween(nuclei,{0,1},0.25)}
                            })};

    const MolecularGeometry A = {
            nuclei,
            ElectronsVector({
                                    {Spin::beta, nuclei.positionsVector()[1]+Eigen::Vector3d(0,0.01,0.03)},
                                    {Spin::beta,inbetween(nuclei,{0,2},0.26)},
                                    {Spin::alpha, nuclei.positionsVector()[2]},
                            })};

    Eigen::VectorXd ekin = Eigen::VectorXd::Zero(3);
    std::vector<Sample> samples = {{A.electrons(),ekin},{B.electrons(),ekin}};


    ParticleKit::create(A);

    General::settings.mode = General::Mode::alchemical;
    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    SOAPClusterer::settings.similarityThreshold = 0.99;
    SOAPClusterer::settings.distanceMatrixCovarianceTolerance = 0.1;
    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;
    Cutoff::settings.radius = 4.0;
    Cutoff::settings.width = 3.0;

    Group gA({1.0, A.electrons(), 0});
    Group gB({2.0, B.electrons(), 1});
    Group maxima({gA,gB});
    maxima.sortAll();

    SOAPClusterer clusterer(nuclei, samples);

    std::cout  << std::endl;
    std::cout << maxima << std::endl;

    clusterer.cluster(maxima);

    std::cout << maxima << std::endl;

    ASSERT_EQ(maxima.size(), 1);
    ASSERT_EQ(maxima[0].size(), 2);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0, 1));
    //ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    //ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1, 2));
    //ASSERT_THAT(maxima[0][1][0].allSampleIds(), ElementsAre(1));
    //ASSERT_THAT(maxima[0][1][1].allSampleIds(), ElementsAre(2));
    //ASSERT_THAT(maxima[0][2].allSampleIds(), ElementsAre(3));


    // Group expected({{A,{B,C},D},{E,F}});
}

class H4test : public ::testing::Test {
public:
    Group A, B, C;
    std::vector<Sample> samples;
    Eigen::VectorXd ekin;

    AtomsVector atoms;

    void SetUp() override {
        //spdlog::set_level(spdlog::level::off);

        atoms = TestMolecules::H4::linear::nuclei.atoms();
        A = Group({1.0, TestMolecules::H4::linear::ionicA.electrons(), 0});
        B = Group({1.0, TestMolecules::H4::linear::ionicAreflected.electrons(), 1});
        C = Group({1.0, TestMolecules::H4::linear::ionicAreflectedReorderedNumbering.electrons(), 2});

        ekin = Eigen::VectorXd::Zero(4);

        samples = {
                {A.representative()->maximum(), ekin},
                {B.representative()->maximum(), ekin},
                {C.representative()->maximum(), ekin}};
    }
};

TEST_F(H4test, TwoClusters_chemical) {
    ParticleKit::create(atoms, A.representative()->maximum());
    SOAPClusterer clusterer(atoms, samples);

    General::settings.mode = General::Mode::chemical;
    SOAPClusterer::settings.similarityThreshold = 1.0;
    SOAPClusterer::settings.distanceMatrixCovarianceTolerance = 0.001;

    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;
    Cutoff::settings.radius = 4.0;
    Cutoff::settings.width = 3.0;

    Group maxima({A, B});

    clusterer.cluster(maxima);

    ASSERT_EQ(maxima.size(), 1);
    ASSERT_EQ(maxima[0].size(), 2);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0, 1));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1));

    // Group expected({{A,B}});
}

TEST_F(H4test, TwoClusters_alchemical) {
    ParticleKit::create(atoms, A.representative()->maximum());
    SOAPClusterer clusterer(atoms, samples);

    General::settings.mode = General::Mode::alchemical;
    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    SOAPClusterer::settings.similarityThreshold = 1.0;
    SOAPClusterer::settings.distanceMatrixCovarianceTolerance = 0.001;

    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;
    Cutoff::settings.radius = 4.0;
    Cutoff::settings.width = 3.0;

    Group maxima({A, B});

    clusterer.cluster(maxima);

    ASSERT_EQ(maxima.size(), 1);
    ASSERT_EQ(maxima[0].size(), 2);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0, 1));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(1));

    // Group expected({{A,B}});
}

TEST_F(H4test, TwoClusters_alchemical_spinswap) {
    ParticleKit::create(atoms, A.representative()->maximum());
    SOAPClusterer clusterer(atoms, samples);

    General::settings.mode = General::Mode::alchemical;
    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 1.0;
    SOAPClusterer::settings.similarityThreshold = 1.0;
    SOAPClusterer::settings.distanceMatrixCovarianceTolerance = 0.001;

    Angular::settings.lmax = 3;
    Radial::settings.nmax = 3;
    Cutoff::settings.radius = 4.0;
    Cutoff::settings.width = 3.0;

    Group maxima({A,C});

    clusterer.cluster(maxima);

    ASSERT_EQ(maxima.size(), 1);
    ASSERT_EQ(maxima[0].size(), 2);
    ASSERT_THAT(maxima[0].allSampleIds(), ElementsAre(0, 2));
    ASSERT_THAT(maxima[0][0].allSampleIds(), ElementsAre(0));
    ASSERT_THAT(maxima[0][1].allSampleIds(), ElementsAre(2));

    // Group expected({{A,C}});
}