/* Copyright (C) 2018-2019 Michael Heuer.
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
#include "CoulombPotential.h"
#include "TestMolecules.h"

using namespace testing;

void compareMatrix(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2){
    ASSERT_EQ(mat1.rows(), mat2.rows());
    ASSERT_EQ(mat1.cols(), mat2.cols());

    for (long i = 0; i < mat1.rows(); ++i)
        for (long j = 0; j < mat1.cols(); ++j)
            ASSERT_EQ(mat1(i,j),mat2(i,j));
}

TEST(ACoulombPotentialTest, energies) {
    auto mol = TestMolecules::CoulombPotentialTest::HeH;
    auto ne = mol.electrons().numberOfEntities();
    auto nn = mol.atoms().numberOfEntities();


    auto Ven = CoulombPotential::energies<Spin,Element>(mol.electrons(),mol.atoms());
    auto Vee = CoulombPotential::energies<Spin>(mol.electrons());
    auto Vnn = CoulombPotential::energies<Element>(mol.atoms());

    auto inf = std::numeric_limits<double>::infinity();

    Eigen::MatrixXd expectedVen(ne,nn);
    expectedVen <<
                -inf,-0.5,\
                -2,-1,\
                -inf,-0.5;
    compareMatrix(Ven,expectedVen);

    Eigen::MatrixXd expectedVee(ne,ne);
    expectedVee <<
                0,1,inf,\
                1,0,1,\
                inf,1,0;
    compareMatrix(Vee,expectedVee);

    Eigen::MatrixXd expectedVnn(nn,nn);
    expectedVnn <<
                0,1,\
                1,0;
    compareMatrix(Vnn,expectedVnn);
}
