//
// Created by leonard on 10.05.19.
//

#include <gmock/gmock.h>
#include <TestMolecules.h>
#include <NearestElectrons.h>
#include <BestMatchDistance.h>

using namespace testing;

class ANearestElectronsTest : public Test {
public:
    const MolecularGeometry &BH3 = TestMolecules::BH3::ionic;
    const ElectronsVector &electrons = BH3.electrons();
    const ElectronsVector &electrons2 = TestMolecules::BH3::ionicMirrored.electrons();
    const ElectronsVector &electrons3 = TestMolecules::BH3::ionicMirrored2.electrons();
    const AtomsVector &nuclei = BH3.atoms();
};

TEST_F(ANearestElectronsTest, GetNonValenceIndices) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electrons, nuclei, 0);
    std::list<long> reference = std::list<long>({0, 4});
    ASSERT_EQ(reference, indices);

    indices = NearestElectrons::getNonValenceIndices(electrons, nuclei, 1);
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
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(electrons, nuclei, position, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByAtomIndex) {
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(electrons, nuclei, 2, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByAtomIndices) {
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(electrons, nuclei, 0, 2, 4);
    std::list<long> reference({0, 2, 4, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByPosition) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electrons, nuclei,
                                                                         nuclei[0].position(), 2);
    std::list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByAtomIndex) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electrons, nuclei, 0, 2);
    std::list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByAtomIndices) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electrons, nuclei, 0, 2, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, PickElements) {
    ElectronsVector reference;
    reference.append(electrons[2]);
    reference.append(electrons[6]);

    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electrons, nuclei, 0, 2, 2);

    ASSERT_EQ(reference, electrons[indices]);
};

TEST_F(ANearestElectronsTest, PickElementsMethod) {
    ElectronsVector reference;
    reference.append(electrons[2]);
    reference.append(electrons[6]);

    ElectronsVector slicedVector = NearestElectrons::getNearestValenceElectrons(electrons, nuclei, 0, 2, 2);
    ASSERT_EQ(reference, slicedVector);
};

TEST_F(ANearestElectronsTest, BestMatch) {
    std::list<long> indices1 = NearestElectrons::getNearestValenceIndices(electrons.positionsVector(), nuclei, 0, 2, 2);
    std::list<long> indices2 = NearestElectrons::getNearestValenceIndices(electrons3.positionsVector(), nuclei, 0, 2, 2);

    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
            electrons[indices1].positionsVector(),
            electrons3[indices2].positionsVector());

    ASSERT_EQ(norm,0);
};

TEST_F(ANearestElectronsTest, BestMatch2) {
    std::list<long> indices1 = NearestElectrons::getNearestValenceIndices(electrons.positionsVector(), nuclei, 0, 1, 2);
    std::list<long> indices2 = NearestElectrons::getNearestValenceIndices(electrons2.positionsVector(), nuclei, 0, 1, 2);

    auto[norm, perm] = BestMatch::Distance::compare<Eigen::Infinity, 2>(
            electrons[indices1].positionsVector(),
            electrons2[indices2].positionsVector());
    ASSERT_NE(norm,0);
};
