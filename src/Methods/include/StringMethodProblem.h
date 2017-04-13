//
// Created by heuer on 13.04.17.
//

#ifndef AMOLQCGUI_STRINGMETHODPROBLEM_H
#define AMOLQCGUI_STRINGMETHODPROBLEM_H

#include "problem.h"
#include "ElectronicWaveFunction.h"
#include "BSpline.h"
#include "ArcLengthParametrizedBSpline.h"

class StringMethodProblem : public cppoptlib::Problem<double,Eigen::Dynamic>{
public:
    StringMethodProblem(long numberOfStates, long numberOfCoords);
    //StringMethodProblem(const Eigen::MatrixXd &initalChainGuess);

    double value(const Eigen::VectorXd &x);

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad);

    bool callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x);

    void resetString();
    void reparametrizeString();
    void constructChainOfStates();

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
    Eigen::MatrixXd chain_;
    Eigen::VectorXd chainGradient_;
    long numberOfStates_, numberOfCoords_;
    unsigned valueCallCount_, gradientCallCount_;
    ElectronicWaveFunction wf_;
    //BSplines::ArcLengthParametrizedBSpline arcLengthParametrizedBSpline_; // use arclength parametrized spline
    BSplines::BSpline bSpline_;
};

#endif //AMOLQCGUI_STRINGMETHODPROBLEM_H
