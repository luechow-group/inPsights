//
// Created by Michael Heuer on 17.08.18.
//

#include <gtest/gtest.h>
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
