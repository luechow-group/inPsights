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

#ifndef INPSIGHTS_LOCALNEWTONSEARCH_H
#define INPSIGHTS_LOCALNEWTONSEARCH_H

#include <solver/newtonraphsonsolver.h>
#include <ElectronicWaveFunctionProblem.h>

#include <meta.h>

class LocalNewtonSearch{
public:

    //LocalNewtonSearch(){};

    template <typename ProblemType>
    void search(ProblemType problem, Eigen::VectorXd &x) const {

        cppoptlib::NewtonRaphsonSolver<ProblemType> solver;
        solver.setDebug(cppoptlib::DebugLevel::High);
        auto crit(cppoptlib::Criteria<double>::defaults());
        crit.gradNorm = 1e-6;
        crit.iterations = 10000;
        solver.setStopCriteria(crit);

        solver.minimize(problem,x);
    };

};

#endif //INPSIGHTS_LOCALNEWTONSEARCH_H
