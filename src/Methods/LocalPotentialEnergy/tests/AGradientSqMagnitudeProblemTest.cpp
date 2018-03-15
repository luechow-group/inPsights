//
// Created by Leonard Reuter on 12.03.18.
//

#include <gtest/gtest.h>
#include <Eigen/Core>
#include "GradientSqMagnitudeProblem.h"
#include "TestProblems.h"
#include "solver/gradientdescentsimplesolver.h"

using namespace testing;
using namespace Eigen;
using namespace TestProblems;

class AGradientSqMagnitudeProblemTest : public Test {};

TEST_F(AGradientSqMagnitudeProblemTest, empty) {
}

TEST_F(AGradientSqMagnitudeProblemTest, Construction) {
    testConstraint problem;
    GradientSqMagnitudeProblem<testConstraint> gradientSqMagnitudeProblem(problem);
}

TEST_F(AGradientSqMagnitudeProblemTest, value) {
    testConstraint problem;

    GradientSqMagnitudeProblem<testConstraint> gradientSqMagnitudeProblem(problem);

    Eigen::VectorXd z(2);
    z << -2,3;

    ASSERT_DOUBLE_EQ(gradientSqMagnitudeProblem.value(z),52);
}

TEST_F(AGradientSqMagnitudeProblemTest, gradient) {
    testConstraint problem;

    GradientSqMagnitudeProblem<testConstraint> gradientSqMagnitudeProblem(problem);

    Eigen::VectorXd z(2);
    z << -2,3;

    Eigen::VectorXd grad = z;
    gradientSqMagnitudeProblem.gradient(z, grad);

    Eigen::VectorXd ref(2);
    ref << -16,24;

    ASSERT_EQ(grad,ref);
}

TEST_F(AGradientSqMagnitudeProblemTest, GradientDescentSimple) {
    testConstraint problem;

    GradientSqMagnitudeProblem<testConstraint> gradientSqMagnitudeProblem(problem);

    Eigen::VectorXd z(2);
    z << -2,3;

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
    crit.gradNorm = 1e-3;

    cppoptlib::GradientDescentSimpleSolver<GradientSqMagnitudeProblem<testConstraint>> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);

    solver.minimize(gradientSqMagnitudeProblem, z);

    Eigen::VectorXd ref(2);
    ref << 0,0;

    ASSERT_GT(1e-3,(z-ref).norm());
}
