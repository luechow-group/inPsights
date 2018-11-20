//
// Created by heuer on 06.04.17.
//

#ifndef INPSIGHTS_STRINGMETHOD_H
#define INPSIGHTS_STRINGMETHOD_H

#include "ElectronicWaveFunction.h"
#include "BSpline.h"
#include "ArcLengthParametrizedBSpline.h"
#include "meta.h"
#include "ChainOfStates.h"
//#include "problemobserver.h"

class StringMethod {//: public cppoptlib::ProblemObserver{

public:
    StringMethod(ChainOfStates chain);
    void optimizeString();

    ChainOfStates & getChain(){ return chain_; }

    BSplines::ArcLengthParametrizedBSpline& getArcLengthParametrizedBSpline(){ return arcLengthParametrizedBSpline_; }

private:
    void minimizeOrthogonalToString();
    void performStep();
    void reparametrizeString(); // specify u values
    void distributeStates();
    void discretizeStringToChain();
    void calculateUnitTangents();

    ElectronicWaveFunction& wf_;
    BSplines::ArcLengthParametrizedBSpline arcLengthParametrizedBSpline_; // use arclength parametrized spline
    ChainOfStates chain_;
    Eigen::MatrixXd unitTangents_;
    Eigen::VectorXd uValues_;
    cppoptlib::Status status_;
};


#endif //INPSIGHTS_STRINGMETHOD_H
