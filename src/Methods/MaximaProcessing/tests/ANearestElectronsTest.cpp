//
// Created by leonard on 10.05.19.
//

#include <gmock/gmock.h>
#include <TestMolecules.h>
#include <NearestElectrons.h>

using namespace testing;
using namespace std;
using namespace NearestElectrons;

class ANearestElectronsTest : public Test {
};

TEST_F(ANearestElectronsTest, GetNonValenceIndices) {
    list<long> indices = getNonValenceIndices(TestMolecules::BH3::ionic, 0);
    list<long> reference = list<long>({0, 4});
    ASSERT_EQ(reference, indices);

    indices = getNonValenceIndices(TestMolecules::BH3::ionic, 1);
    reference = list<long>({});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetNonValenceIndicesAll) {
    list<long> indices = getNonValenceIndices(TestMolecules::BH3::ionic);
    list<long> reference({0, 4});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByPosition) {
    Eigen::Vector3d position = TestMolecules::inbetween(TestMolecules::BH3::ionic.atoms(), {0, 2}, 0.7);
    list<long> indices = getNearestElectronsIndices(TestMolecules::BH3::ionic, position, 2);
    list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByAtomIndex) {
    list<long> indices = getNearestElectronsIndices(TestMolecules::BH3::ionic, 2, 2);
    list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByAtomIndices) {
    list<long> indices = getNearestElectronsIndices(TestMolecules::BH3::ionic, 0, 2, 4);
    list<long> reference({0, 2, 4, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByPosition) {
    list<long> indices = getNearestValenceIndices(TestMolecules::BH3::ionic,
                                                  TestMolecules::BH3::ionic.atoms()[0].position(), 2);
    list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByAtomIndex) {
    list<long> indices = getNearestValenceIndices(TestMolecules::BH3::ionic, 0, 2);
    list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByAtomIndices) {
    list<long> indices = getNearestValenceIndices(TestMolecules::BH3::ionic, 0, 2, 2);
    list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};
