//
// Created by heuer on 19.10.18.
//

#include <gmock/gmock.h>
#include <BestMatchDistanceDensityBasedClusterer.h>
#include <BestMatchDistanceSimilarityClusterer.h>
#include <random>
#include <TestMolecules.h>
#include <Sample.h>
#include <Reference.h>
#include <NaturalConstants.h>
#include <Metrics.h>

using namespace testing;

class ABestMatchDistanceDensityBasedClustererTest : public ::testing::Test {
public:
    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        BestMatchDistanceSimilarityClusterer::settings.similarityRadius = 0.1; // prevent assert
    }
    Group makeRingLikeCluster(Group &references, std::vector<Sample> &samples, unsigned n){
        auto rng = std::default_random_engine(static_cast<unsigned long>(std::clock()));
        const auto &normal = TestMolecules::fourElectrons::normal.electrons();

        // emplace first element (the first elements determines the ordering of the cluster)
        references.emplace_back(Reference(0,normal));

        Eigen::PermutationMatrix<Eigen::Dynamic> perm(normal.numberOfEntities());
        perm.setIdentity();

        // start with second element
        bool anyWrong = false;
        while (!anyWrong){
            for (unsigned i = 1; i < n; ++i) {
                double angle = Constant::pi/2.0*double(i)/double(n); // 90°
                Eigen::Vector3d zAxis = {0,0,1};
                auto evCopy = normal;
                evCopy.positionsVector().rotateAroundOrigin(angle, zAxis);

                // random permutation
                std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), rng);
                evCopy.permute(perm);
                references.emplace_back(Reference(std::pow(std::sin(angle),2),evCopy));

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

TEST_F(ABestMatchDistanceDensityBasedClustererTest, RotationallySymmetricCluster){
    BestMatchDistanceDensityBasedClusterer::settings.clusterRadius = 0.1;

    unsigned n = 20;
    Group references;
    std::vector<Sample> samples;
    makeRingLikeCluster(references, samples, n);


    auto rng = std::default_random_engine(static_cast<unsigned long>(std::clock()));

    ASSERT_EQ(references.size(), n);
    std::shuffle(references.begin(), references.end(), rng);

    BestMatchDistanceDensityBasedClusterer globalClusterSorter(samples);
    globalClusterSorter.cluster(references);

    ASSERT_EQ(references.size(),1);
    ASSERT_EQ(references[0].size(), n);

    auto referenceDistanceMatrix = Metrics::positionalDistances(references[0][0].representative()->maximum().positionsVector());

    for (size_t i = 1; i < references[0].size(); ++i) {
        auto distanceMatrix = Metrics::positionalDistances(references[0][i].representative()->maximum().positionsVector());
        ASSERT_TRUE(distanceMatrix.isApprox(referenceDistanceMatrix));
    }
}


TEST_F(ABestMatchDistanceDensityBasedClustererTest, RotationallySymmetricAndPointLikeCluster){
    BestMatchDistanceDensityBasedClusterer::settings.clusterRadius = 0.1;

    unsigned n = 20;
    Group references;
    std::vector<Sample> samples;
    makeRingLikeCluster(references, samples, n);

    const auto &ionic = TestMolecules::fourElectrons::ionic.electrons();

    Eigen::PermutationMatrix<Eigen::Dynamic> perm(ionic.numberOfEntities());
    perm.setIdentity();

    auto rng = std::default_random_engine(static_cast<unsigned long>(std::clock()));

    references.emplace_back(Reference(10, ionic));
    unsigned m = 5;
    for (unsigned i = 1; i < m; ++i) {
        auto maxDev = 1./std::sqrt(3.0)*BestMatchDistanceDensityBasedClusterer::settings.clusterRadius.get();

        Eigen::Vector3d randomTranslation;
        std::uniform_real_distribution<double> uniformRealDistribution(-maxDev, maxDev);
        randomTranslation.x() = uniformRealDistribution(rng);
        randomTranslation.y() = uniformRealDistribution(rng);
        randomTranslation.z() = uniformRealDistribution(rng);

        auto evCopy = ionic;
        evCopy.positionsVector().translate(randomTranslation);

        // random permutation
        std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), rng);
        evCopy.permute(perm);
        references.emplace_back(Reference(10 ,evCopy));

        Sample s(evCopy, Eigen::VectorXd::Random(ionic.numberOfEntities()));
        samples.emplace_back(std::move(s));
    }

    ASSERT_EQ(references.size(), n+m);

    std::shuffle(references.begin(), references.end(), rng);

    BestMatchDistanceDensityBasedClusterer globalClusterSorter(samples);
    globalClusterSorter.cluster(references);

    ASSERT_EQ(references.size(), 2);
    ASSERT_EQ(references[0].size(), n);
    ASSERT_EQ(references[1].size(), m);

    for (size_t j = 0; j < references.size(); j++){
        auto referenceDistanceMatrix = Metrics::positionalDistances(references[j][0].representative()->maximum().positionsVector());

        for (size_t i = 1; i < references[j].size(); ++i) {
            auto distanceMatrix = Metrics::positionalDistances(references[j][i].representative()->maximum().positionsVector());
            ASSERT_TRUE(distanceMatrix.isApprox(referenceDistanceMatrix));
        }
    }
}
