//
// Created by heuer on 06.04.17.
//

#include "StringMethod.h"

template<typename ProblemType, typename SolverType, int Ord>
void StringMethod<ProblemType,SolverType,Ord>::resetString(unsigned numberOfStates) {
    //intialize solvers
    solvers_.clear();

    for (unsigned i = 0; i < numberOfStates; ++i) {
        SolverType solver;
        solvers_.push_back(solver);
    }

}

template<typename ProblemType,typename SolverType, int Ord>
void StringMethod<ProblemType,SolverType,Ord>::evaluateString(Eigen::VectorXd &x) {
    double f = problemReference.value(x);
    std::cout << f << std::endl;
}
