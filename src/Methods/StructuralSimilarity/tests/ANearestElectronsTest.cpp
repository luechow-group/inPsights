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
    const PositionsVector &electronPositions = BH3.electrons().positionsVector();
    const AtomsVector &nuclei = BH3.atoms();
};

TEST_F(ANearestElectronsTest, GetNonValenceIndices) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electronPositions, nuclei, 0);
    std::list<long> reference = std::list<long>({0, 4});
    ASSERT_EQ(reference, indices);

    indices = NearestElectrons::getNonValenceIndices(electronPositions, nuclei, 1);
    reference = std::list<long>({});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetNonValenceIndicesAll) {
    std::list<long> indices = NearestElectrons::getNonValenceIndices(electronPositions, nuclei);
    std::list<long> reference({0, 4});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByPosition) {
    Eigen::Vector3d position = TestMolecules::inbetween(nuclei, {0, 2}, 0.7);
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(electronPositions, nuclei, position, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByAtomIndex) {
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(electronPositions, nuclei, 2, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetElectronsByAtomIndices) {
    std::list<long> indices = NearestElectrons::getNearestElectronsIndices(electronPositions, nuclei, 0, 2, 4);
    std::list<long> reference({0, 2, 4, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByPosition) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electronPositions, nuclei,
                                                                         nuclei[0].position(), 2);
    std::list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByAtomIndex) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electronPositions, nuclei, 0, 2);
    std::list<long> reference({6, 7});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, GetValenceByAtomIndices) {
    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electronPositions, nuclei, 0, 2, 2);
    std::list<long> reference({2, 6});
    ASSERT_EQ(reference, indices);
};

TEST_F(ANearestElectronsTest, ConvertIndexList) {
    Eigen::ArrayXi reference(6);
    reference << 6, 7, 8, 18, 19, 20;

    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electronPositions, nuclei, 0, 2, 2);

    Eigen::ArrayXi indicesArray = NearestElectrons::longIndexListToPositionArrayXi(indices);
    ASSERT_EQ(reference.matrix(), indicesArray.matrix());
};

TEST_F(ANearestElectronsTest, SliceVector) {
    Eigen::VectorXd reference(6);
    reference << electronPositions[2], electronPositions[6];

    std::list<long> indices = NearestElectrons::getNearestValenceIndices(electronPositions, nuclei, 0, 2, 2);
    Eigen::ArrayXi indicesArray = NearestElectrons::longIndexListToPositionArrayXi(indices);

    // indexArray.unaryExpr(x) does the same as python/fortran x(indexArray)
    Eigen::VectorXd positionsVector = indicesArray.unaryExpr(electronPositions.asEigenVector());

    ASSERT_EQ(reference, positionsVector);
};

TEST_F(ANearestElectronsTest, SliceVectorMethod) {
    Eigen::VectorXd reference(6);
    reference << electronPositions[2], electronPositions[6];

    Eigen::VectorXd positionsVector = NearestElectrons::getNearestValencePositions(electronPositions, nuclei, 0, 2, 2);
    ASSERT_EQ(reference, positionsVector);
};
