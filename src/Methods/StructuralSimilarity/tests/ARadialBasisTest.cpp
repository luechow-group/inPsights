//
// Created by dahl on 21.11.17.
//

#include <gtest/gtest.h>
#include <Eigen/Core>
#include <RadialBasis.h>
#include <iomanip>

using namespace testing;

class ARadialBasisTest : public ::testing::Test {};

TEST_F(ARadialBasisTest, CorrectInverseMatrixSqrt) {

    int nmax = 15;
    double rcutoff = 2.0;

    RadialBasis radialBasis(nmax, rcutoff);

    Eigen::Matrix2d mat;
    mat << 4,1,1,4;

    Eigen::Matrix2d reference;
    reference <<
              +(1.0/(2.0*std::sqrt(3.0))) + 1.0/(2.0*std::sqrt(5.0)),\
              -(1.0/(2.0*std::sqrt(3.0))) + 1.0/(2.0*std::sqrt(5.0)),\
              -(1.0/(2.0*std::sqrt(3.0))) + 1.0/(2.0*std::sqrt(5.0)),\
              +(1.0/(2.0*std::sqrt(3.0))) + 1.0/(2.0*std::sqrt(5.0));

    Eigen::Matrix2d invSqrtMat = radialBasis.inverseMatrixSqrt(mat);



    ASSERT_TRUE((reference * reference).inverse().isApprox(mat));
    ASSERT_TRUE(invSqrtMat.isApprox(reference));
    ASSERT_TRUE((invSqrtMat * invSqrtMat).inverse().isApprox(mat));

}


