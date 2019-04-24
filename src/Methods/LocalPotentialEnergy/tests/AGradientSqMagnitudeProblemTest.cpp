/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <gmock/gmock.h>
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
    TestConstraint problem;
    GradientSqMagnitudeProblem<TestConstraint> gradientSqMagnitudeProblem(problem);
}

TEST_F(AGradientSqMagnitudeProblemTest, value) {
    TestConstraint problem;

    GradientSqMagnitudeProblem<TestConstraint> gradientSqMagnitudeProblem(problem);

    Eigen::VectorXd z(2);
    z << -2,3;

    ASSERT_DOUBLE_EQ(gradientSqMagnitudeProblem.value(z),52);
}

TEST_F(AGradientSqMagnitudeProblemTest, gradient) {
    TestConstraint problem;

    GradientSqMagnitudeProblem<TestConstraint> gradientSqMagnitudeProblem(problem);

    Eigen::VectorXd z(2);
    z << -2,3;

    Eigen::VectorXd grad = z;
    gradientSqMagnitudeProblem.gradient(z, grad);

    Eigen::VectorXd ref(2);
    ref << -16,24;

    ASSERT_EQ(grad,ref);
}

TEST_F(AGradientSqMagnitudeProblemTest, GradientDescentSimple) {
    TestConstraint problem;

    GradientSqMagnitudeProblem<TestConstraint> gradientSqMagnitudeProblem(problem);

    Eigen::VectorXd z(2);
    z << -2,3;

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
    crit.gradNorm = 1e-3;

    cppoptlib::GradientDescentSimpleSolver<GradientSqMagnitudeProblem<TestConstraint>> solver;
    solver.setDebug(cppoptlib::DebugLevel::None);
    solver.setStopCriteria(crit);

    solver.minimize(gradientSqMagnitudeProblem, z);

    Eigen::VectorXd ref(2);
    ref << 0,0;

    ASSERT_GT(1e-3,(z-ref).norm());
}
