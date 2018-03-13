//
// Created by Leonard Reuter on 12.03.18.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include "solver/gradientdescentsolver.h"
#include "TestProblems.h"

using namespace testing;
using namespace Eigen;
using namespace TestProblems;

class ATestProblemsTest : public Test {};

TEST_F(ATestProblemsTest, testProblemValue) {
    testProblem problem;

    Eigen::VectorXd z(2);
    z << -2,3;

    ASSERT_DOUBLE_EQ(problem.value(z),1);
}

TEST_F(ATestProblemsTest, testConstraintValue) {
    testConstraint constraint;

    Eigen::VectorXd z(2);
    z << -2,3;

    ASSERT_DOUBLE_EQ(constraint.value(z),13);
}

TEST_F(ATestProblemsTest, testProblemGradient) {
    testProblem problem;

    Eigen::VectorXd z(2);
    z << -2,3;

    Eigen::VectorXd grad = z;

    problem.gradient(z,grad);

    Eigen::VectorXd reference(2);
    reference << 1,1;
    ASSERT_EQ(grad,reference);
}

TEST_F(ATestProblemsTest, testConstraintGradient) {
    testConstraint constraint;

    Eigen::VectorXd z(2);
    z << -2,3;

    Eigen::VectorXd grad = z;

    constraint.gradient(z,grad);

    Eigen::VectorXd reference(2);
    reference << -4,6;

    ASSERT_DOUBLE_EQ(grad[0],reference[0]);
    ASSERT_DOUBLE_EQ(grad[1],reference[1]);
}