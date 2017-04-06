//
// Created by heuer on 06.04.17.
//

#ifndef AMOLQCGUI_STRINGMETHOD_H
#define AMOLQCGUI_STRINGMETHOD_H

#include <Eigen/Core>
#include <vector>
#include "ElectronicWaveFunctionProblem.h"

class StringMethod {
    //using Scalar = typename ProblemType_::Scalar;
    //using ProblemType = ProblemType_;

public:
    StringMethod( ElectronicWaveFunctionProblem &objFunc )
            : objFunc_( objFunc ) {}

    void evaluateString(Eigen::VectorXd &x);

private:
    ElectronicWaveFunctionProblem &objFunc_;
};

#endif //AMOLQCGUI_STRINGMETHOD_H
