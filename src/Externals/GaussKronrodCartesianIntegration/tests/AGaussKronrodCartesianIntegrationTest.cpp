//
// Created by Michael Heuer on 12.03.18.
//

#include <gtest/gtest.h>

#include "GaussKronrodCartesianIntegration.h"
#include "IntegrandExampleFunction.h"

using namespace testing;

class AGaussKronrodCartesianIntegrationTest : public ::testing::Test {};

TEST_F(AGaussKronrodCartesianIntegrationTest, NumericalIntegrationExample) {

    double alpha = 1;
    IntegrandExampleFunctor f(alpha);

    GaussKronrod::Integrator<double> integrator(500);

    GaussKronrod::Integrator<double>::QuadratureRule quadratureRule = Eigen::Integrator<double>::GaussKronrod15;

    // Define the desired absolute and relative errors.
    double desAbsErr = 0;
    double desRelErr = Eigen::NumTraits<double>::epsilon() * 50;

    IntegrandExampleFunctor * fptr= &f;

    double result = integrator.quadratureAdaptive(*fptr, 0, 1, desAbsErr, desRelErr, quadratureRule);
    double expected = -4.0;


    double absErr = std::abs(expected-result);
    double relErr = absErr/expected;

    ASSERT_NEAR(result, expected, absErr);
    ASSERT_LE(relErr, desRelErr);
}

