/* Copyright (C) 2018-2019 Michael Heuer.
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
#include <Reference.h>
#include <Sample.h>
#include <TestMolecules.h>
#include <IdentityClusterer.h>
#include <PreClusterer.h>
#include <algorithm>
#include <random>

using namespace testing;

class AIdentityClustererTest : public ::testing::Test {
public:
    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        PreClusterer::settings.radius = 10; // prevent assert
    }

    std::pair<Cluster, std::vector<Sample>> prepareRandomTripletMaxima(std::default_random_engine &rng) {
        Cluster tripletMaxima;
        std::vector<Sample> tripletSamples;
        AtomsVector atoms = TestMolecules::H2::ElectronsInCores::normal.atoms();

        // Multiplicity = 3: Spin flip is not possible.
        tripletMaxima = {
                Cluster({atoms, 1.00, ElectronsVector({{Spin::alpha, {0, 0, 1.00}},
                                              {Spin::alpha, {0, 0, 0}}}), 0}),
                Cluster({atoms, 1.01, ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                              {Spin::alpha, {0, 0, 1.01}}}), 1}),
                Cluster({atoms, 1.02, ElectronsVector({{Spin::alpha, {0, 0, 1.02}},
                                              {Spin::alpha, {0, 0, 0}}}), 2}),
                Cluster({atoms, 1.03, ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                              {Spin::alpha, {0, 0, 1.03}}}), 3}),
                Cluster({atoms, 1.10, ElectronsVector({{Spin::alpha, {0, 0, 1.10}},
                                              {Spin::alpha, {0, 0, 0}}}), 4}),
                Cluster({atoms, 1.11, ElectronsVector({{Spin::alpha, {0, 0, 1.11}},
                                              {Spin::alpha, {0, 0, 0}}}), 5}),
                Cluster({atoms, 1.12, ElectronsVector({{Spin::alpha, {0, 0, 1.12}},
                                              {Spin::alpha, {0, 0, 0}}}), 6}),
                Cluster({atoms, 1.13, ElectronsVector({{Spin::alpha, {0, 0, 1.13}},
                                              {Spin::alpha, {0, 0, 0}}}), 7}),
        };
        std::shuffle(std::begin(tripletMaxima), std::end(tripletMaxima), rng);
        for (auto &i : tripletMaxima) {
            Sample s(ElectronsVector({{i.representative()->maximum().typesVector()[0], {0, 0, 0}},
                                      {i.representative()->maximum().typesVector()[1], {0, 0, 0}}}),
                     Eigen::VectorXd::Random(2));
            tripletSamples.emplace_back(std::move(s));
        }
        tripletMaxima.sortAll();
        return {tripletMaxima, tripletSamples};
    }

    std::pair<Cluster, std::vector<Sample>> prepareRandomSingletMaxima(std::default_random_engine &rng) {
        Cluster singletMaxima;
        std::vector<Sample> singletSamples;
        AtomsVector atoms = TestMolecules::H2::ElectronsInCores::normal.atoms();

        // Multiplicity = 1: Spin flip is possible.
        singletMaxima = {
                Cluster({atoms, 1.00, ElectronsVector({{Spin::alpha, {0, 0, 1.00}},
                                              {Spin::beta,  {0, 0, 0}}}), 0}),
                Cluster({atoms, 1.01, ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                              {Spin::beta,  {0, 0, 1.01}}}), 1}),
                Cluster({atoms, 1.02, ElectronsVector({{Spin::alpha, {0, 0, 1.02}},
                                              {Spin::beta,  {0, 0, 0}}}), 2}),
                Cluster({atoms, 1.03, ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                              {Spin::beta,  {0, 0, 1.03}}}), 3}),
                Cluster({atoms, 1.10, ElectronsVector({{Spin::alpha, {0, 0, 1.10}},
                                              {Spin::beta,  {0, 0, 0}}}), 4}),
                Cluster({atoms, 1.11, ElectronsVector({{Spin::alpha, {0, 0, 1.11}},
                                              {Spin::beta,  {0, 0, 0}}}), 5}),
                Cluster({atoms, 1.12, ElectronsVector({{Spin::alpha, {0, 0, 1.12}},
                                              {Spin::beta,  {0, 0, 0}}}), 6}),
                Cluster({atoms, 1.13, ElectronsVector({{Spin::alpha, {0, 0, 1.13}},
                                              {Spin::beta,  {0, 0, 0}}}), 7}),
        };
        std::shuffle(std::begin(singletMaxima), std::end(singletMaxima), rng);


        for (auto &i : singletMaxima) {
            Sample s(ElectronsVector({{i.representative()->maximum().typesVector()[0], {0, 0, 0}},
                                      {i.representative()->maximum().typesVector()[1], {0, 0, 0}}}),
                     Eigen::VectorXd::Random(2));
            singletSamples.emplace_back(std::move(s));
        }
        singletMaxima.sortAll();
        return {singletMaxima, singletSamples};
    }
};

TEST_F(AIdentityClustererTest, OneListTriplet) {
    IdentityClusterer::settings.radius = 2;
    IdentityClusterer::settings.valueIncrement = 1;

    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for (auto seed : std::vector<unsigned long>{0, randomSeed}) {
        auto rng = std::default_random_engine(seed);

        auto triplet = prepareRandomTripletMaxima(rng);
        IdentityClusterer globalIdentiySorter(triplet.second);
        globalIdentiySorter.cluster(triplet.first);

        ASSERT_THAT(triplet.first.representative()->sampleIds(), ElementsAre(0, 1, 2, 3, 4, 5, 6, 7));
    }
}

TEST_F(AIdentityClustererTest, OneListSinglet) {
    IdentityClusterer::settings.radius = 2;
    IdentityClusterer::settings.valueIncrement = 1;

    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for (auto seed : std::vector<unsigned long>{0, randomSeed}) {
        auto rng = std::default_random_engine(seed);

        auto singlet = prepareRandomSingletMaxima(rng);
        IdentityClusterer globalIdentiySorter(singlet.second);
        globalIdentiySorter.cluster(singlet.first);

        ASSERT_THAT(singlet.first.representative()->sampleIds(), ElementsAre(0, 1, 2, 3, 4, 5, 6, 7));
    }
}

TEST_F(AIdentityClustererTest, TwoListsTriplet) {
    IdentityClusterer::settings.radius = 1;
    IdentityClusterer::settings.valueIncrement = 0.05;

    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for (auto seed : std::vector<unsigned long>{0, randomSeed}) {
        auto rng = std::default_random_engine(seed);

        auto triplet = prepareRandomTripletMaxima(rng);
        IdentityClusterer globalIdentiySorter(triplet.second);
        globalIdentiySorter.cluster(triplet.first);

        ASSERT_THAT(triplet.first.at(0).representative()->sampleIds(), ElementsAre(0, 1, 2, 3));
        ASSERT_THAT(triplet.first.at(1).representative()->sampleIds(), ElementsAre(4, 5, 6, 7));
    }
}

TEST_F(AIdentityClustererTest, TwoListsSinglet) {
    IdentityClusterer::settings.radius = 1;
    IdentityClusterer::settings.valueIncrement = 0.05;

    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for (auto seed : std::vector<unsigned long>{0, randomSeed}) {
        auto rng = std::default_random_engine(seed);

        auto singlet = prepareRandomSingletMaxima(rng);
        IdentityClusterer globalIdentiySorter(singlet.second);
        globalIdentiySorter.cluster(singlet.first);

        ASSERT_THAT(singlet.first.at(0).representative()->sampleIds(), ElementsAre(0, 1, 2, 3));
        ASSERT_THAT(singlet.first.at(1).representative()->sampleIds(), ElementsAre(4, 5, 6, 7));
    }
}

TEST_F(AIdentityClustererTest, TwoListsIncrementBorderCaseTriplet) {
    IdentityClusterer::settings.radius = 1;
    IdentityClusterer::settings.valueIncrement = 0.1;

    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for (auto seed : std::vector<unsigned long>{0, randomSeed}) {
        auto rng = std::default_random_engine(seed);

        auto triplet = prepareRandomTripletMaxima(rng);
        IdentityClusterer globalIdentiySorter(triplet.second);
        globalIdentiySorter.cluster(triplet.first);

        ASSERT_THAT(triplet.first.at(0).representative()->sampleIds(), ElementsAre(0, 1, 2, 3, 4));
        ASSERT_THAT(triplet.first.at(1).representative()->sampleIds(), ElementsAre(5, 6, 7));
    }
}

TEST_F(AIdentityClustererTest, TwoListsIncrementBorderCaseSinglet) {
    IdentityClusterer::settings.radius = 1;
    IdentityClusterer::settings.valueIncrement = 0.1;

    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for (auto seed : std::vector<unsigned long>{0, randomSeed}) {
        auto rng = std::default_random_engine(seed);

        auto singlet = prepareRandomSingletMaxima(rng);
        IdentityClusterer globalIdentiySorter(singlet.second);
        globalIdentiySorter.cluster(singlet.first);

        ASSERT_THAT(singlet.first.at(0).representative()->sampleIds(), ElementsAre(0, 1, 2, 3, 4));
        ASSERT_THAT(singlet.first.at(1).representative()->sampleIds(), ElementsAre(5, 6, 7));
    }
}
