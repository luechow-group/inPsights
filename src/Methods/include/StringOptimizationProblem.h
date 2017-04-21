//
// Created by heuer on 13.04.17.
//

#ifndef AMOLQCGUI_STRINGMETHODPROBLEM_H
#define AMOLQCGUI_STRINGMETHODPROBLEM_H

#include "problem.h"
#include "StringMethod.h"

class StringOptimizationProblem : public cppoptlib::Problem<double,Eigen::Dynamic>{
public:

    /*
     StringMethod(long numberOfStates, long numberOfCoords,
     ElectronicWaveFunction & wf,
     Eigen::VectorXd unitTangent
     )
     */

    StringOptimizationProblem(long numberOfStates,
                              long numberOfCoords,
                              ElectronicWaveFunction &wf,
                              const Eigen::VectorXd &unitTangent
                              //Eigen::MatrixXd &chain
                              );
    //StringMethodProblem(const Eigen::MatrixXd &initalChainGuess);

    double value(const Eigen::VectorXd &x);

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &orthogonalGrad);

    bool callback(cppoptlib::Criteria<double> &state, Eigen::VectorXd &x);


    //void setChainOfStates(const Eigen::VectorXd &systemCoordVector);

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
    unsigned stepCounter_;
    long numberOfStates_, numberOfCoords_;
    //Eigen::MatrixXd &chain_;
    ElectronicWaveFunction &wf_;
    const Eigen::VectorXd &unitTangent_;
    unsigned valueCallCount_, gradientCallCount_;


};

#endif //AMOLQCGUI_STRINGMETHODPROBLEM_H
