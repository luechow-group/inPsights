//
// Created by heuer on 06.04.17.
//

#ifndef AMOLQCGUI_STRINGMETHOD_H
#define AMOLQCGUI_STRINGMETHOD_H

#include <Eigen/Core>
#include <vector>
#include "ElectronicWaveFunctionProblem.h"
#include "ProblemObserver.h"

class StringMethod : public ProblemObserver{
    //using Scalar = typename ProblemType_::Scalar;
    //using ProblemType = ProblemType_;

public:
    StringMethod( ElectronicWaveFunctionProblem &problem )
            : problemReference( problem ) {}

    void resetString();

    void evaluateString(Eigen::VectorXd &x);

    void stepPerformed() override;


private:
    ElectronicWaveFunctionProblem &problemReference;
};

#endif //AMOLQCGUI_STRINGMETHOD_H
