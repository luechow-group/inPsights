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

using namespace testing;

class ABestMatchDistanceDensityBasedClustererTest : public ::testing::Test {
public:
    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
        BestMatchDistanceSimilarityClusterer::settings.similarityRadius = 0.1; // prevent assert
    }
};

TEST_F(ABestMatchDistanceDensityBasedClustererTest, RotationallySymmetricCluster){
    BestMatchDistanceDensityBasedClusterer::settings.clusterRadius = 0.1;

    const auto &normal = TestMolecules::threeElectrons::normal.electrons();

    std::vector<Sample> samples;
    Group references, expected;

    auto rng = std::default_random_engine(static_cast<unsigned long>(std::clock()));

    // emplace first element (the first elements determines the ordering of the cluster)
    expected.emplace_back(Reference(0,normal));
    references.emplace_back(Reference(0,normal));

    Eigen::PermutationMatrix<Eigen::Dynamic> perm(normal.numberOfEntities());
    perm.setIdentity();

    unsigned n = 20;
    // start with second element
    for (unsigned i = 1; i < n; ++i) {
        double angle = Constant::pi*double(i)/double(n); // 180°
        Eigen::Vector3d xAxis = {1,0,0};
        auto evCopy = normal;
        evCopy.positionsVector().rotateAroundOrigin(angle, xAxis);
        expected.emplace_back(Reference(std::pow(std::sin(angle),2),evCopy));

        // random permuation
        std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), rng);
        evCopy.permute(perm);
        references.emplace_back(Reference(std::pow(std::sin(angle),2),evCopy));

        Sample s(evCopy, Eigen::VectorXd::Random(normal.numberOfEntities()));
        samples.emplace_back(std::move(s));
    }

    ASSERT_EQ(references.size(), expected.size());

    // don't shuffle the first element
    std::shuffle(references.begin()+1, references.end(), rng);

    BestMatchDistanceDensityBasedClusterer globalClusterSorter(samples);
    globalClusterSorter.cluster(references);

    ASSERT_EQ(references.size(),1);
    ASSERT_EQ(references[0].size(), expected.size());

    for (size_t i = 0; i < references[0].size(); ++i) {
        ASSERT_TRUE(references[0][i].representative()->maximum() == expected[i].representative()->maximum());
        ASSERT_TRUE(references[0][i].representative()->value() == expected[i].representative()->value());
    }
}


TEST_F(ABestMatchDistanceDensityBasedClustererTest, RotationallyAndPointSymmetricCluster){
    BestMatchDistanceDensityBasedClusterer::settings.clusterRadius = 0.1;

    const auto &normal = TestMolecules::threeElectrons::normal.electrons();
    auto ionic = TestMolecules::threeElectrons::normal.electrons();
    ionic.positionsVector().translate({1.0,0,0});

    assert(normal.numberOfEntities() == ionic.numberOfEntities());

    std::vector<Sample> samples;
    Group references, cluster0Expected, cluster1Expected;

    auto rng = std::default_random_engine(static_cast<unsigned long>(std::clock()));

    // emplace first element
    cluster0Expected.emplace_back(Reference(0,normal));
    references.emplace_back(Reference(0,normal));


    Eigen::PermutationMatrix<Eigen::Dynamic> perm(normal.numberOfEntities());
    perm.setIdentity();

    unsigned n = 20;
    // start with second element
    for (unsigned i = 1; i < n; ++i) {
        double angle = Constant::pi*double(i)/double(n); // 180°
        Eigen::Vector3d xAxis = {1,0,0};
        auto evCopy = normal;
        evCopy.positionsVector().rotateAroundOrigin(angle, xAxis);
        cluster0Expected.emplace_back(Reference(std::pow(std::sin(angle),2),evCopy));

        // random permuation
        std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), rng);
        evCopy.permute(perm);
        references.emplace_back(Reference(std::pow(std::sin(angle),2),evCopy));

        Sample s(evCopy, Eigen::VectorXd::Random(normal.numberOfEntities()));
        samples.emplace_back(std::move(s));
    }
    ASSERT_EQ(references.size(), n);
    ASSERT_EQ(references.size(), cluster0Expected.size());

    references.emplace_back(Reference(0, ionic));
    cluster1Expected.emplace_back(Reference(0, ionic));

    unsigned m = 5;
    for (unsigned i = 1; i < m; ++i) {
        auto maxDev = 1./std::sqrt(3.0)*BestMatchDistanceSimilarityClusterer::settings.similarityRadius.get();

        Eigen::Vector3d randomTranslation; //= Eigen::Vector3d::Random(-maxDev, maxDev);
        std::uniform_real_distribution<double> uniformRealDistribution(-maxDev, maxDev);
        randomTranslation.x() = uniformRealDistribution(rng);
        randomTranslation.y() = uniformRealDistribution(rng);
        randomTranslation.z() = uniformRealDistribution(rng);

        auto evCopy = ionic;
        evCopy.positionsVector().translate(randomTranslation);
        cluster0Expected.emplace_back(Reference(0, evCopy));

        // random permuation
        std::shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size(), rng);
        evCopy.permute(perm);
        references.emplace_back(Reference(0 ,evCopy));

        Sample s(evCopy, Eigen::VectorXd::Random(normal.numberOfEntities()));
        samples.emplace_back(std::move(s));
    }

    ASSERT_EQ(references.size(), n+m);
    ASSERT_EQ(references.size(), cluster0Expected.size() + cluster1Expected.size());
    ASSERT_EQ(references[0].size(), cluster0Expected[0].size());
    ASSERT_EQ(references[1].size(), cluster0Expected[1].size());

    // don't shuffle the first element
    std::shuffle(references.begin()+1, references.end(), rng);

    BestMatchDistanceDensityBasedClusterer globalClusterSorter(samples);
    globalClusterSorter.cluster(references);

    for (size_t i = 0; i < references[0].size(); ++i) {
        ASSERT_TRUE(references[0][i].representative()->maximum() == cluster0Expected[i].representative()->maximum());
        ASSERT_TRUE(references[0][i].representative()->value() == cluster0Expected[i].representative()->value());
    }
}
