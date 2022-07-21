// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <DensityBasedClusterer.h>
#include <PreClusterer.h>
#include <random>
#include <TestMolecules.h>
#include <Sample.h>
#include <Maximum.h>
#include <NaturalConstants.h>
#include <Metrics.h>

using namespace testing;

class ADensityBasedClustererTest : public ::testing::Test {
public:
    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        PreClusterer::settings.radius = 0.1; // prevent assert
    }

    Cluster makeRingLikeCluster(Cluster &references, std::vector<Sample> &samples, unsigned n, std::default_random_engine& rng){
        const auto &normal = TestMolecules::fourElectrons::normal;

        // emplace first element (the first elements determines the ordering of the cluster)
        references.emplace_back(Maximum(normal.atoms(), 0, normal.electrons(), 0));

        // start with second element
        bool anyWrong = false;
        while (!anyWrong){
            for (unsigned i = 1; i < n; ++i) {
                double angle = Constant::pi/2.0*double(i)/double(n); // 90°
                Eigen::Vector3d zAxis = {0,0,1};
                auto evCopy = normal.electrons();
                evCopy.positionsVector().rotateAroundOrigin(angle, zAxis);

                // random permutation
                evCopy.permute(evCopy.randomPermutation(rng));
                references.emplace_back(Maximum(normal.atoms(), std::pow(std::sin(angle), 2), evCopy, 0));

                Sample s(evCopy, Eigen::VectorXd::Random(normal.numberOfEntities()));
                samples.emplace_back(std::move(s));
            }

            auto referenceDistanceMatrix2 = Metrics::positionalDistances(references[0].representative()->maximum().positionsVector());

            // test if any permutation is wrong
            for (size_t i = 1; i < references.size(); ++i) {
                auto distanceMatrix2 = Metrics::positionalDistances(references[i].representative()->maximum().positionsVector());
                if (!distanceMatrix2.isApprox(referenceDistanceMatrix2)){
                    anyWrong = true;
                    break;
                };
            }
        }
        return references;
    }
};

TEST_F(ADensityBasedClustererTest, RotationallySymmetricCluster){
    DensityBasedClusterer::settings.radius = 0.1;

    unsigned n = 20;

    // check different seeds
    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for(auto seed : std::vector<unsigned long>{0,randomSeed}) {
        auto rng = std::default_random_engine(seed);

        Cluster references;
        std::vector<Sample> samples;

        makeRingLikeCluster(references, samples, n, rng);

        ASSERT_EQ(references.size(), n);
        std::shuffle(references.begin(), references.end(), rng);

        DensityBasedClusterer globalClusterSorter(samples);
        globalClusterSorter.cluster(references);

        ASSERT_EQ(references.size(), 1);
        ASSERT_EQ(references[0].size(), n);

        auto referenceDistanceMatrix = Metrics::positionalDistances(
                references[0][0].representative()->maximum().positionsVector());

        for (size_t i = 1; i < references[0].size(); ++i) {
            auto distanceMatrix = Metrics::positionalDistances(
                    references[0][i].representative()->maximum().positionsVector());
            ASSERT_TRUE(distanceMatrix.isApprox(referenceDistanceMatrix));
        }
    }
}

TEST_F(ADensityBasedClustererTest, RotationallySymmetricAndPointLikeCluster){
    DensityBasedClusterer::settings.radius = 0.1;

    unsigned n = 20;
    unsigned m = 5;

    const auto &ionic = TestMolecules::fourElectrons::ionic;

    // check different seeds
    auto randomSeed = static_cast<unsigned long>(std::clock());
    std::cout << "random seed: " << randomSeed << std::endl;

    for(auto seed : std::vector<unsigned long>{0,randomSeed}) {
        auto rng = std::default_random_engine(seed);

        Cluster references;
        std::vector<Sample> samples;
        makeRingLikeCluster(references, samples, n, rng);

        references.emplace_back(Maximum(ionic.atoms(), 10, ionic.electrons(), 0));
        for (unsigned i = 1; i < m; ++i) {
            auto evCopy = ionic.electrons();
            evCopy.positionsVector().shake(DensityBasedClusterer::settings.radius.get(), rng);

            // random permutation
            evCopy.permute(evCopy.randomPermutation(rng));
            references.emplace_back(Maximum(ionic.atoms(), 10, evCopy, 0));

            Sample s(evCopy, Eigen::VectorXd::Random(ionic.numberOfEntities()));
            samples.emplace_back(std::move(s));
        }

        ASSERT_EQ(references.size(), n + m);

        std::shuffle(references.begin(), references.end(), rng);

        DensityBasedClusterer globalClusterSorter(samples);
        globalClusterSorter.cluster(references);

        ASSERT_EQ(references.size(), 2);
        ASSERT_EQ(references[0].size(), n);
        ASSERT_EQ(references[1].size(), m);

        for (size_t j = 0; j < references.size(); j++) {
            auto referenceDistanceMatrix =
                    Metrics::positionalDistances(references[j][0].representative()->maximum().positionsVector());

            for (size_t i = 1; i < references[j].size(); ++i) {
                auto distanceMatrix =
                        Metrics::positionalDistances(references[j][i].representative()->maximum().positionsVector());
                // for the comparison, the threshold has to be 2*radius
                // since every particle can be displaced by radius, the distance between
                // two particles can change by 2*radius
                ASSERT_TRUE(distanceMatrix.isApprox(referenceDistanceMatrix,
                                                    2 * DensityBasedClusterer::settings.radius.get()));
            }
        }
    }
}
