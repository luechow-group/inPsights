//
// Created by Leonard Reuter on 12.03.18.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include "LagrangeProblem.h"
#include "solver/gradientdescentsolver.h"

using namespace testing;
using namespace Eigen;

class ALagrangeProblemTest : public Test {public:
    void SetUp() override {};
};

class testProblem: public cppoptlib::Problem<double,Eigen::Dynamic>{
public:
    double value(const Eigen::VectorXd & x) override {
        double value = 0.0;
        for (int i = 0; i < x.size(); i++){
            value += x[i];
        }
    };

    void gradient(const Eigen::VectorXd & x, Eigen::VectorXd &grad) override {
        for (int i = 0; i < x.size(); i++){
            grad[i] = 1;
        }
    };
};

class testConstraint : public cppoptlib::Problem<double,Eigen::Dynamic>{
public:
    double value(const Eigen::VectorXd & x) override {
        double value = 0.0;
        for (int i = 0; i < x.size(); i++){
            value += x[i]*x[i];
        }
    };

    void gradient(const Eigen::VectorXd & x, Eigen::VectorXd &grad) override {
        for (int i = 0; i < x.size(); i++){
            grad[i] = 2*x[i];
        }
    };
};

TEST(ALagrangeProblemTest, Construction) {
    testProblem problem;
    testConstraint constraint;
    LagrangeProblem<testProblem,testConstraint> lagrangeProblem(problem,constraint,1);
}

//took test functions from wikipedia
TEST(ALagrangeProblemTest, Value) {
    testProblem problem;
    testConstraint constraint;
    LagrangeProblem<testProblem,testConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd z(2);
    z << -2,-2;

    ASSERT_DOUBLE_EQ(lagrangeProblem.getProblem().value(z),-4);
    ASSERT_DOUBLE_EQ(lagrangeProblem.getConstraint().value(z),8);
    ASSERT_DOUBLE_EQ(lagrangeProblem.value(y),10);
}

TEST(ALagrangeProblemTest, Gradient) {
    testProblem problem;
    testConstraint constraint;
    LagrangeProblem<testProblem,testConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd gradient = y;
    lagrangeProblem.gradient(y,gradient);

    Eigen::VectorXd reference(3);
    reference << -7,-7,7;
    ASSERT_EQ(gradient, reference);
}

TEST(ALagrangeProblemTest, Optimization) {
    testProblem problem;
    testConstraint constraint;
    LagrangeProblem<testProblem,testConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd gradient = y;
    lagrangeProblem.gradient(y,gradient);

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();

    cppoptlib::GradientDescentSolver<LagrangeProblem<testProblem,testConstraint>> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);

    solver.minimize(lagrangeProblem, y);

    Eigen::VectorXd reference(3);
    reference << -0.70713277,-0.70713277,0.70713277;

    ASSERT_GT(1e-8,(y-reference).norm());
}

TEST(ALagrangeProblemTest, getProblem) {
    testProblem problem;
    testConstraint constraint;
    LagrangeProblem<testProblem,testConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd z(2);
    z << -2,-2;

    ASSERT_EQ(lagrangeProblem.getProblem().value(z), problem.value(z));
}

TEST(ALagrangeProblemTest, getConstraint) {
    testProblem problem;
    testConstraint constraint;
    LagrangeProblem<testProblem,testConstraint> lagrangeProblem(problem,constraint,1);

    Eigen::VectorXd y(3);
    y << -2,-2,2;

    Eigen::VectorXd z(2);
    z << -2,-2;

    ASSERT_EQ(lagrangeProblem.getConstraint().value(z), constraint.value(z));
}