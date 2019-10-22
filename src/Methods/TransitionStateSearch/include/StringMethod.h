/* Copyright (C) 2017-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
