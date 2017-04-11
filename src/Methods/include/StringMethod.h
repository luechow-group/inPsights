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

    StringMethod( ProblemType &problem, unsigned numberOfStates)
            : numberOfStates(numberOfStates),
              problemReference(problem) {}

    void resetString(unsigned numberOfStates);

    void evaluateString(Eigen::VectorXd &x);

    void stepPerformed() override{
        std::cout << "step performed" << std::endl;
    };

private:
    unsigned numberOfStates;
    ProblemType& problemReference;
    std::vector<cppoptlib::ISolver<ProblemType,Ord>> solvers_;
};

#endif //AMOLQCGUI_STRINGMETHOD_H
