//
// Created by heuer on 06.04.17.
//

#ifndef AMOLQCGUI_STRINGMETHOD_H
#define AMOLQCGUI_STRINGMETHOD_H

#include <Eigen/Core>
#include <vector>
#include "problem.h"
#include "problemobserver.h"
#include "solver/isolver.h"

template<typename ProblemType, typename SolverType, int Ord>
class StringMethod : public cppoptlib::ProblemObserver{
public:
    using Scalar    = typename ProblemType::Scalar;
    using TVector   = typename ProblemType::TVector;

    StringMethod(unsigned numberOfStates)
            : numberOfStates(numberOfStates) {
        problem.addObserver(this);
        resetString(numberOfStates);
    }

    void resetString(unsigned numberOfStates) {
        solvers_.clear();
        for (unsigned i = 0; i < numberOfStates; ++i) {
            solvers_.push_back(SolverType());
        }
    }

    void evaluateString(Eigen::VectorXd &x) {
        double f = problem.value(x);
        std::cout << f << std::endl;
    }

    void stepPerformed() override{
        std::cout << "step performed" << std::endl;
    };

private:
    unsigned numberOfStates;
    ProblemType problem;
    std::vector<SolverType> solvers_;
};

#endif //AMOLQCGUI_STRINGMETHOD_H
