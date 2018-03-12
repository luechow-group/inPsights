//
// Created by Leonard Reuter on 12.03.18.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include <iostream>
#include "PotentialProblem.h"
#include "ElectronsVector.h"


using namespace testing;
using namespace Eigen;

class APotentialProblemTest : public Test {public:
    void SetUp() override {
        Atom atom1(Vector3d(1,2,3),Elements::ElementType::Ag);
        Atom atom2(Vector3d(-1,0,3.5),Elements::ElementType::Au);
        Atom atom3(Vector3d(-7.3,0.5,9),Elements::ElementType::C);

        atomsVector.append(atom1);
        atomsVector.append(atom2);
        atomsVector.append(atom3);

        Electron elec1(Vector3d(1,2,3.1));
        Electron elec2(Vector3d(-1,0,3.6));
        Electron elec3(Vector3d(-7.3,0.5,9.1));

        electronsVector.append(elec1);
        electronsVector.append(elec2);
        electronsVector.append(elec3);

        atomsVector2.append(Atom(Vector3d(0,0,0),Elements::ElementType::H));
        electronsVector2.append(Electron(Vector3d(0,0,2)));
    }
    AtomsVector atomsVector, atomsVector2;
    ElectronsVector electronsVector, electronsVector2;
};

TEST_F(APotentialProblemTest, DefaultConstruction) {
    PotentialProblem potentialProblem();
}

TEST_F(APotentialProblemTest, Construction) {
    PotentialProblem potentialProblem(atomsVector);
}

TEST_F(APotentialProblemTest, getAtoms) {
    PotentialProblem potentialProblem(atomsVector);

    ASSERT_EQ(atomsVector.positionsVector().positionsAsEigenVector(),
              potentialProblem.getAtoms().positionsVector().positionsAsEigenVector());
    ASSERT_EQ(atomsVector.elementTypesVector().elementTypesAsEigenVector(),
              potentialProblem.getAtoms().elementTypesVector().elementTypesAsEigenVector());
}

TEST_F(APotentialProblemTest, Value1e1c) {
    PotentialProblem potentialProblem(atomsVector2);
    VectorXd x = electronsVector2.positionsVector().positionsAsEigenVector();

    ASSERT_DOUBLE_EQ(potentialProblem.value(x),-0.5);
}

TEST_F(APotentialProblemTest, Gradient1e1c) {
    PotentialProblem potentialProblem(atomsVector2);
    VectorXd x = electronsVector2.positionsVector().positionsAsEigenVector();

    auto grad = x;
    potentialProblem.gradient(x,grad);
    ASSERT_DOUBLE_EQ(grad[2],0.25);
}

TEST_F(APotentialProblemTest, Value3e3c) {
    PotentialProblem potentialProblem(atomsVector);

    VectorXd x = electronsVector.positionsVector().positionsAsEigenVector();
    ASSERT_DOUBLE_EQ(potentialProblem.value(x),-2.0005632341974531);
}