//
// Created by heuer on 13.04.17.
//

#ifndef AMOLQCGUI_STRINGMETHODPROBLEM_H
#define AMOLQCGUI_STRINGMETHODPROBLEM_H

#include "problem.h"
#include "ElectronicWaveFunction.h"

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

};

#endif //AMOLQCGUI_STRINGMETHODPROBLEM_H
