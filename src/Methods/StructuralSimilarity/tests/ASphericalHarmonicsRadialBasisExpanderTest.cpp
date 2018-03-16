//
// Created by Michael Heuer on 16.03.18.
//

#include <gtest/gtest.h>
#include <LebedevSphericalIntegration/TestFunctions.h>
#include <LebedevSphericalIntegration/GridCreator.h>
#include "SphericalHarmonicsRadialBasisExpander.h"

#include "iostream"
using namespace testing;

class ASphericalHarmonicsRadialBasisExpanderTest : public ::testing::Test {};

TEST_F(ASphericalHarmonicsRadialBasisExpanderTest, CoefficientsVector) {
    unsigned nmax = 15; //TODO 20 does not work, why? => Matrix to big for inversion
    unsigned lmax = 2;
    double rCutoff = 5.0;

    Gauss3d f;
    SphericalHarmonicsRadialBasisExpander expander(nmax,lmax,rCutoff);

    auto fExp = ExpandedFunction(f,nmax,lmax,rCutoff);

    Eigen::Vector3d r({1,2,3});

    Lebedev::GridCreator gridCreator;
    gridCreator.changeGrid(Lebedev::OrderType::LD0014);

    Eigen::MatrixX3d grid = gridCreator.grid().leftCols(3);
    std::cout << grid << std::endl;

    //TODO also iterate over r
    for (int i = 0; i < grid.rows(); ++i) {
        ASSERT_NEAR(f(grid.row(i)), fExp(grid.row(i)),1e-4);
    }


}