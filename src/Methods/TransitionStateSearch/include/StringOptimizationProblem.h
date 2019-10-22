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

#ifndef INPSIGHTS_STRINGMETHODPROBLEM_H
#define INPSIGHTS_STRINGMETHODPROBLEM_H

#include "problem.h"
#include "ElectronicWaveFunction.h"
#include <fstream>

enum StateGradientType{ SimpleGradient, OrthogonalToString, Fixed, ClimbingImage};

class StringOptimizationProblem : public cppoptlib::Problem<double,Eigen::Dynamic>{
public:


    StringOptimizationProblem(long numberOfStates,
                              long numberOfCoords,
                              ElectronicWaveFunction &wf,
                              Eigen::MatrixXd &unitTangents //TODO make const
                              //Eigen::MatrixXd &chain
                              );
    //StringMethodProblem(const Eigen::MatrixXd &initalChainGuess);

    double value(const Eigen::VectorXd &x) override;

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) override;

    bool callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) override;

    Eigen::VectorXd stateValues(const Eigen::VectorXd &x);

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
    
    void fixGradient(Eigen::VectorXd& grad);

    Eigen::VectorXd getNucleiPositions() const;

    void putElectronsIntoNuclei(Eigen::VectorXd& x, Eigen::VectorXd& grad);

    std::vector<unsigned long> getIndicesOfElectronsNotAtNuclei();
    std::vector<unsigned long> getIndicesOfElectronsAtNuclei();

private:
    unsigned stepCounter_;
    long numberOfStates_, numberOfCoords_;
    ElectronicWaveFunction &wf_;
    Eigen::MatrixXd &unitTangents_; //TODO make const
    Eigen::VectorXd values_;
    unsigned valueCallCount_, gradientCallCount_;

    std::vector<StateGradientType> stateTypes_;

    std::vector<unsigned long> indicesOfElectronsNotAtNuclei_;
    std::vector<unsigned long> indicesOfElectronsAtNuclei_;
    std::ofstream outfile;
};

#endif //INPSIGHTS_STRINGMETHODPROBLEM_H
