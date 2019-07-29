//
// Created by leonard on 10.05.19.
//

#include <gmock/gmock.h>
#include <TestMolecules.h>
#include <NearestElectrons.h>
#include <Metrics.h>

using namespace testing;

class ANearestElectronsTest : public Test {
public:
    const MolecularGeometry &BH3 = TestMolecules::BH3::ionic;
    const ElectronsVector &electrons = BH3.electrons();
    const ElectronsVector &electrons2 = TestMolecules::BH3::ionicRotated.electrons();
    ElectronsVector electrons3;
    const AtomsVector &nuclei = BH3.atoms();

    void SetUp() override {
        electrons3 =
                ElectronsVector({
                    electrons2[0],
                    electrons2[1],
                    electrons2[2],
                    electrons2[7],
                    electrons2[4],
                    electrons2[5],
                    electrons2[6],
                    electrons2[3]
                });
    }
    std::list<long>
    getNearestElectronsIndices(const ElectronsVector &electrons, const AtomsVector &nuclei,
                             const Eigen::Vector3d &position,
                             const long &count){
        std::function<double(const Eigen::Vector3d &, const std::vector<Eigen::Vector3d> &)> distanceFunction = Metrics::minimalDistance<2>;
        return NearestElectrons::getNearestElectronsIndices(electrons, nuclei, std::vector<Eigen::Vector3d>({position}), count, 100.0,
                                          distanceFunction, true);
    };
};

TEST_F(ANearestElectronsTest, GetNonValenceIndices) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electrons, nuclei[0]);
    std::list<long> reference = std::list<long>({0, 4});
    ASSERT_EQ(reference, indices);

    indices = NearestElectrons::getNonValenceIndices(electrons, nuclei[1]);
    reference = std::list<long>({});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetNonValenceIndicesAll) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electrons, nuclei);
    std::list<long> reference({0, 4});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByPosition) {
    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.7);
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(electrons, position, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByPosition) {
    std::list<long> indices = getNearestElectronsIndices(electrons, nuclei,
                                                                         nuclei[0].position(), 2);
    std::list<long> reference({7, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, PickElements) {
    ElectronsVector reference;
    reference.append(electrons[6]);
    reference.append(electrons[2]);

    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.5);

    std::list<long> indices = getNearestElectronsIndices(electrons, nuclei, position, 2);

    ASSERT_EQ(reference, electrons[indices]);
};