//
// Created by Leonard Reuter on 12.03.18.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include "LagrangeProblem.h"
#include "solver/gradientdescentsolver.h"
#include "solver/gradientdescentsimplesolver.h"
#include "solver/bfgssolver.h"
#include "TestProblems.h"

using namespace testing;
using namespace Eigen;
using namespace TestProblems;

class ALagrangeProblemTest : public Test {};

TEST_F(ALagrangeProblemTest, Construction) {
    TestProblem problem;
    TestConstraint constraint;
    LagrangeProblem<TestProblem,TestConstraint> lagrangeProblem(problem,constraint,1);
}

//took test functions from wikipedia
TEST_F(ALagrangeProblemTest, Value) {
    TestProblem problem;
    TestConstraint constraint;
    LagrangeProblem<TestProblem,TestConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd z(2);
    z << -2,-2;

    ASSERT_DOUBLE_EQ(lagrangeProblem.getProblem().value(z),-4);
    ASSERT_DOUBLE_EQ(lagrangeProblem.getConstraint().value(z),8);
    ASSERT_DOUBLE_EQ(lagrangeProblem.value(y),10);
}

TEST_F(ALagrangeProblemTest, Gradient) {
    TestProblem problem;
    TestConstraint constraint;
    LagrangeProblem<TestProblem,TestConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd gradient = y;
    lagrangeProblem.gradient(y,gradient);

    Eigen::VectorXd reference(3);
    reference << -7,-7,7;
    ASSERT_EQ(gradient, reference);
}

TEST_F(ALagrangeProblemTest, GradientDescentSolver) {
    TestProblem problem;
    TestConstraint constraint;
    LagrangeProblem<TestProblem,TestConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd gradient = y;
    lagrangeProblem.gradient(y,gradient);

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();

    cppoptlib::GradientDescentSolver<LagrangeProblem<TestProblem,TestConstraint>> solver;
    solver.setDebug(cppoptlib::DebugLevel::None);
    solver.setStopCriteria(crit);

    solver.minimize(lagrangeProblem, y);

    Eigen::VectorXd reference(3);
    reference << -0.70713277,-0.70713277,0.70713277;

    ASSERT_GT(1e-8,(y-reference).norm());
}

TEST_F(ALagrangeProblemTest, GradientDescentSimpleSolver) {
    TestProblem problem;
    TestConstraint constraint;
    LagrangeProblem<TestProblem,TestConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd gradient = y;
    lagrangeProblem.gradient(y,gradient);

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
    crit.gradNorm = 1e-3;

    cppoptlib::GradientDescentSimpleSolver<LagrangeProblem<TestProblem,TestConstraint>> solver;
    solver.setDebug(cppoptlib::DebugLevel::None);
    solver.setStopCriteria(crit);

    solver.minimize(lagrangeProblem, y);

    Eigen::VectorXd reference(3);
    reference << -0.70713277,-0.70713277,0.70713277;

    ASSERT_GT(1e-3,(y-reference).norm());
}


TEST_F(ALagrangeProblemTest, BfgsSolver) {
    TestProblem problem;
    TestConstraint constraint;
    LagrangeProblem<TestProblem,TestConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -0.7,-0.7,0.7;

    Eigen::VectorXd gradient = y;
    lagrangeProblem.gradient(y,gradient);

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
    crit.gradNorm = 1e-9;

    cppoptlib::BfgsSolver<LagrangeProblem<TestProblem,TestConstraint>> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);

    solver.minimize(lagrangeProblem, y);

    Eigen::VectorXd reference(3);
    reference << -0.70713277,-0.70713277,0.70713277;

    ASSERT_GT(1e-4,(y-reference).norm());
}

TEST_F(ALagrangeProblemTest, getProblem) {
    TestProblem problem;
    TestConstraint constraint;
    LagrangeProblem<TestProblem,TestConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd z(2);
    z << -2,-2;

    ASSERT_EQ(lagrangeProblem.getProblem().value(z), problem.value(z));
}

TEST_F(ALagrangeProblemTest, getConstraint) {
    TestProblem problem;
    TestConstraint constraint;
    LagrangeProblem<TestProblem,TestConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd z(2);
    z << -2,-2;

    ASSERT_EQ(lagrangeProblem.getConstraint().value(z), constraint.value(z));
}