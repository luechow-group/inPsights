//
// Created by Leonard Reuter on 12.03.18.
//

#include <gtest/gtest.h>
#include <Eigen/Core>
#include "ElectronicWaveFunctionProblem.h"

using namespace testing;
using namespace Eigen;

class AElectronicWaveFunctionProblemTest : public Test {public:
    void SetUp() override {
        /*Atom atom1(Vector3d(1,2,3),Element::Ag);
        Atom atom2(Vector3d(-1,0,3.5),Element::Au);
        Atom atom3(Vector3d(-7.3,0.5,9),Element::C);

        atomsVector.append(atom1);
        atomsVector.append(atom2);
        atomsVector.append(atom3);

        Electron elec1(Vector3d(1,2,3.1));
        Electron elec2(Vector3d(-1,0,3.6));
        Electron elec3(Vector3d(-7.3,0.5,9.1));

        electronsVector.append(elec1);
        electronsVector.append(elec2);
        electronsVector.append(elec3);

        atomsVector2.append(Atom(Vector3d(0,0,0),Element::H));
        electronsVector2.append(Electron(Vector3d(0,0,2)));*/
    }
    /*
    AtomsVector atomsVector, atomsVector2;
    ElectronsVector electronsVector, electronsVector2;*/
};

TEST_F(AElectronicWaveFunctionProblemTest, DefaultConstruction) {
    ElectronicWaveFunctionProblem electronicWaveFunctionProblem;

    std::stringstream oss;
    oss << electronicWaveFunctionProblem.getAtomsVector();
    ASSERT_EQ(oss.str(),"");
}

