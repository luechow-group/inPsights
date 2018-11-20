//
// Created by Michael Heuer on 06.04.18.
//

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
