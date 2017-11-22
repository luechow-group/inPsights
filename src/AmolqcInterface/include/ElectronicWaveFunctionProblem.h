//
// Created by heuer on 06.04.17.
//

#ifndef AMOLQCGUI_ELECTRONICWAVEFUNCTIONPROBLEM_H
#define AMOLQCGUI_ELECTRONICWAVEFUNCTIONPROBLEM_H

#include "ElectronicWaveFunction.h"
#include "problem.h"
#include "ElectronCollections.h"

class ElectronicWaveFunctionProblem : public cppoptlib::Problem<double,Eigen::Dynamic>
{
public:
    explicit ElectronicWaveFunctionProblem(const std::string &fileName);

    double value(const Eigen::VectorXd &x) override;

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) override;

    void hessian(const Eigen::VectorXd&x, Eigen::MatrixXd &hessian) override;

    bool callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x) override;

    Eigen::VectorXd getNucleiPositions();

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

    ElectronCollections getOptimizationPath(){
        return optimizationPath_;
    }

private:
    unsigned valueCallCount_, gradientCallCount_;
    ElectronicWaveFunction& wf_;
    ElectronCollections optimizationPath_;
};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTIONPROBLEM_H
