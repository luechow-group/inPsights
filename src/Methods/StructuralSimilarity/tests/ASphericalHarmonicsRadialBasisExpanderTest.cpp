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
    unsigned nmax = 20;
    unsigned lmax = 0;
    double rCutoff = 5.0;

    Gauss3d f;
    SphericalHarmonicsRadialBasisExpander expander(nmax,lmax,rCutoff);

    auto fExp = ExpandedFunction(f,nmax,lmax,rCutoff);


    Lebedev::GridCreator gridCreator;
    gridCreator.changeGrid(Lebedev::OrderType::LD0014);
    Eigen::MatrixX3d grid = gridCreator.grid().leftCols(3);

    int nDistances = 5;

    // don't start at zero distance
    for (int j = 1; j <= nDistances; ++j) {

        //TODO also iterate over r
        for (int i = 0; i < grid.rows(); ++i) {
            Eigen::Vector3d r = rCutoff*double(j)/double(nDistances)*grid.row(i);
            std::cout << r.transpose() << " " << f(r) << " " << fExp(r) << std::endl;
            //ASSERT_NEAR(f(r), fExp(r),1e-4);
        }
    }

}