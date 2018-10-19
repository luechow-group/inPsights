//
// Created by heuer on 19.10.18.
//

#include <gmock/gmock.h>
#include <SimilarReferences.h>
#include <GlobalClusterSorter.h>
#include <random>
#include <TestMolecules.h>

using namespace testing;

class AGlobalClusterSorterTest : public ::testing::Test {
public:
};

TEST(AGlobalClusterSorterTest, TwoRotatedAndOneStationaryElectrons){
    const auto &ev = TestMolecules::threeElectrons::normal.electrons();

    std::vector<Sample> samples;
    std::vector<Reference> references, expected;

    auto rng = std::default_random_engine{};
    rng.seed(std::clock());

    // emplace first element
    expected.emplace_back(Reference(0,ev));
    references.emplace_back(Reference(0,ev));

    std::vector<int> permIndices;
    for (int i = 0; i < ev.numberOfEntities(); ++i)
        permIndices.emplace_back(i);

    unsigned n = 20;
    // start with second element
    for (unsigned i = 1; i < n-1; ++i) {
        double angle = 1*M_PI*double(i)/double(n-1); // 180°
        Eigen::Vector3d xAxis = {1,0,0};
        auto evCopy = ev;
        evCopy.positionsVector().rotateAroundOrigin(angle, xAxis);
        expected.emplace_back(Reference(std::pow(std::sin(angle),2),evCopy));

        // random permuation
        std::shuffle(permIndices.begin(),permIndices.end(),rng);
        Eigen::Map<Eigen::VectorXi> indices(permIndices.data(), ev.numberOfEntities());
        Eigen::PermutationMatrix<Eigen::Dynamic> perm(indices);
        evCopy.permute(perm);
        references.emplace_back(Reference(std::pow(std::sin(angle),2),evCopy));

        Sample s(evCopy, Eigen::VectorXd::Random(ev.numberOfEntities()));
        samples.emplace_back(std::move(s));
    }

    // don't shuffle the first element
    std::shuffle(references.begin()+1, references.end(), rng);

    std::vector<SimilarReferences> globallySimilarMaxima({SimilarReferences(references.begin())});
    for(auto it = references.begin()+1; it != references.end(); it++){
        globallySimilarMaxima.emplace_back(SimilarReferences(it));
    }

    std::vector<std::vector<SimilarReferences>> globallyClusteredMaxima;
    GlobalClusterSorter globalClusterSorter(samples, globallySimilarMaxima, globallyClusteredMaxima, 1);
    globalClusterSorter.sort();


    ASSERT_EQ(globallyClusteredMaxima.size(),1);
    ASSERT_EQ(globallyClusteredMaxima[0].size(),expected.size());

    for (size_t i = 0; i < globallyClusteredMaxima[0].size(); ++i) {
        ASSERT_TRUE(globallyClusteredMaxima[0][i].representativeReference().maximum()
        == expected[i].maximum());

        ASSERT_TRUE(globallyClusteredMaxima[0][i].representativeReference().value()
        == expected[i].value());
    }
}
