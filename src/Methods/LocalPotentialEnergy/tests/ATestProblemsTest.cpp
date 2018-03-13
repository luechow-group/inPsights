//
// Created by Leonard Reuter on 12.03.18.
//

#include <gtest/gtest.h>
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

TEST_F(ATestProblemsTest, GradientDescent) {
    testConstraint constraint;

    Eigen::VectorXd z(2);
    z << -2,3;

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();

    cppoptlib::GradientDescentSolver<testConstraint> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);

    solver.minimize(constraint, z);

    Eigen::VectorXd ref(2);
    ref << 0,0;

    ASSERT_GT(1e-4,(z-ref).norm());
}