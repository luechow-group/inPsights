/* Copyright (C) 2018 Leonard Reuter.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <gmock/gmock.h>
#include <Eigen/Core>
#include <iostream>
#include "PotentialProblem.h"
#include <ParticlesVector.h>


using namespace testing;
using namespace Eigen;

class APotentialProblemTest : public Test {public:
    void SetUp() override {
        Atom atom1(Element::Ag,{1,2,3});
        Atom atom2(Element::Au,{-1,0,3.5});
        Atom atom3(Element::C,{-7.3,0.5,9});

        atomsVector.append(atom1);
        atomsVector.append(atom2);
        atomsVector.append(atom3);

        Electron elec1({1,2,3.1});
        Electron elec2({-1,0,3.6});
        Electron elec3({-7.3,0.5,9.1});

        electronsVector.append(elec1);
        electronsVector.append(elec2);
        electronsVector.append(elec3);

        atomsVector2.append(Atom(Element::H,{0,0,0}));
        electronsVector2.append(Electron({0,0,2}));
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

    ASSERT_EQ(atomsVector.positionsVector().asEigenVector(),
              potentialProblem.getAtoms().positionsVector().asEigenVector());
    ASSERT_EQ(atomsVector.typesVector().asEigenVector(),
              potentialProblem.getAtoms().typesVector().asEigenVector());
}

TEST_F(APotentialProblemTest, Value1e1c) {
    PotentialProblem potentialProblem(atomsVector2);
    VectorXd x = electronsVector2.positionsVector().asEigenVector();

    ASSERT_DOUBLE_EQ(potentialProblem.value(x),-0.5);
}

TEST_F(APotentialProblemTest, Gradient1e1c) {
    PotentialProblem potentialProblem(atomsVector2);
    VectorXd x = electronsVector2.positionsVector().asEigenVector();

    auto grad = x;
    potentialProblem.gradient(x,grad);
    ASSERT_DOUBLE_EQ(grad[2],0.25);
}

TEST_F(APotentialProblemTest, Value3e3c) {
    PotentialProblem potentialProblem(atomsVector);

    VectorXd x = electronsVector.positionsVector().asEigenVector();
    ASSERT_DOUBLE_EQ(potentialProblem.value(x),-2.0005632341974531);
}