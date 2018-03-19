//
// Created by dahl on 21.11.17.
//

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <RadialBasis.h>
#include <iomanip>

using namespace testing;

class ARadialBasisTest : public ::testing::Test {};

//TEST_F(ARadialBasisTest, CheckWMatrix) {
//    unsigned nmax = 3;
//    Eigen::MatrixXd reference(nmax,nmax);
//    reference << \
//               16.5329,-29.0783, 13.3085,\
//              -29.0783, 65.2568,-36.0001,\
//               13.3085,-36.0001, 23.4921; //
//
//    RadialBasis radialBasis(nmax, 0);
//    auto calculated = radialBasis.W(nmax);
//
//
//    for (int i = 0; i < nmax; ++i) {
//        for (int j = 0; j < nmax; ++j) {
//            auto relError = std::abs(calculated.row(i)[j]- reference.row(i)[j])/reference.row(i)[j];
//            /*std::cout << std::scientific << std::setprecision(18)
//                      << calculated.row(i)[j] << ", "
//                      << reference.row(i)[j]
//                      << relError << std::endl;*/
//            ASSERT_LE(relError,2e-6);
//        }
//    }
//    ASSERT_TRUE(calculated.isApprox(reference,2e-6));
//}

