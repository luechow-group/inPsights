//
// Created by leonard on 10.05.19.
//

#include <gmock/gmock.h>
#include <TestMolecules.h>
#include <NearestElectrons.h>

using namespace testing;

class ANearestElectronsTest : public Test {
public:
    const MolecularGeometry &BH3 = TestMolecules::BH3::ionic;
};

TEST_F(ANearestElectronsTest, GetNonValenceIndices) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(BH3, 0);
    std::list<long> reference = std::list<long>({0, 4});
    ASSERT_EQ(reference, indices);

    indices = NearestElectrons::getNonValenceIndices(BH3, 1);
    reference = std::list<long>({});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetNonValenceIndicesAll) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(BH3);
    std::list<long> reference({0, 4});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByPosition) {
    Eigen::Vector3d position = TestMolecules::inbetween(BH3.atoms(), {0, 2}, 0.7);
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(BH3, position, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByAtomIndex) {
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(BH3, 2, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByAtomIndices) {
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(BH3, 0, 2, 4);
    std::list<long> reference({0, 2, 4, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByPosition) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(BH3,
                                                                         BH3.atoms()[0].position(), 2);
    std::list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByAtomIndex) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(BH3, 0, 2);
    std::list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByAtomIndices) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(BH3, 0, 2, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, ConvertIndexList) {
    Eigen::ArrayXi reference(6);
    reference << 6, 7, 8, 18, 19, 20;

    std::list<long> indices = NearestElectrons::getNearestValenceIndices(BH3, 0, 2, 2);

    Eigen::ArrayXi indicesArray = NearestElectrons::LongIndexListToPositionArrayXi(indices);
    ASSERT_EQ(reference.matrix(), indicesArray.matrix());
};

TEST_F(ANearestElectronsTest, SliceVector) {
    Eigen::VectorXd reference(6);
    reference << BH3.electrons()[2].position(), BH3.electrons()[6].position();

    std::list<long> indices = NearestElectrons::getNearestValenceIndices(BH3, 0, 2, 2);
    Eigen::ArrayXi indicesArray = NearestElectrons::LongIndexListToPositionArrayXi(indices);

    // indexArray.unaryExpr(x) does the same as python/fortran x(indexArray)
    Eigen::VectorXd positionsVector = indicesArray.unaryExpr(BH3.electrons().positionsVector().asEigenVector());

    ASSERT_EQ(reference, positionsVector);
};
