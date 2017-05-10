//
// Created by heuer on 06.04.17.
//

#ifndef AMOLQCGUI_ELECTRONICWAVEFUNCTIONPROBLEM_H
#define AMOLQCGUI_ELECTRONICWAVEFUNCTIONPROBLEM_H

#include <iomanip>
#include "ElectronicWaveFunction.h"
#include "problem.h"
#include "observableproblem.h"
#include "problemobserver.h"

class ElectronicWaveFunctionProblem : public cppoptlib::Problem<double,Eigen::Dynamic>,
                                      public cppoptlib::ObservableProblem
{
public:

    ElectronicWaveFunctionProblem();

    double value(const Eigen::VectorXd &x);

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad);

    bool callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x);

    unsigned getValueCallCount(){
        return valueCallCount_;
    }

    unsigned getGradientCallCount(){
        return gradientCallCount_;
    }

    unsigned getTotalElocCalls(){
        return getValueCallCount()+getGradientCallCount();
    }

    void resetCounters(){
        valueCallCount_ = 0;
        gradientCallCount_ = 0;
    }

private:
    unsigned valueCallCount_, gradientCallCount_;
    ElectronicWaveFunction wf_;
};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTIONPROBLEM_H
