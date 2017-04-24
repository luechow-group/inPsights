//
// Created by heuer on 13.04.17.
//

#include "StringOptimizationProblem.h"
#include <iostream>
#include <iomanip>



StringOptimizationProblem::StringOptimizationProblem(long numberOfStates,
                                                     long numberOfCoords,
                                                     ElectronicWaveFunction &wf,
                                                     Eigen::MatrixXd &unitTangents //TODO make const
                                                     //Eigen::MatrixXd &chain
                                                     )
        : stepCounter_(0),
          numberOfStates_(numberOfStates),
          numberOfCoords_(numberOfCoords),
          wf_(wf),
          unitTangents_(unitTangents),
          //chain_(chain),
          valueCallCount_(0),
          gradientCallCount_(0)
{
    assert(numberOfStates_ > 2);
    assert( (numberOfCoords_ > 0) && (numberOfCoords_%3 == 0) );
}


double StringOptimizationProblem::value(const Eigen::VectorXd &x) {

    double value = 0;
    for (int i = 0; i < numberOfStates_; ++i) {
        Eigen::VectorXd xi(x.segment(i*numberOfCoords_, numberOfCoords_) );

        valueCallCount_++;
        wf_.evaluate(xi);

        value += wf_.getNegativeLogarithmizedProbabilityDensity();
        //chain_(i,0) = wf_.getNegativeLogarithmizedProbabilityDensity();
    }
    return value;//chain_.col(0).sum();
}

void StringOptimizationProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    Eigen::VectorXd stateGrad(numberOfCoords_);
    Eigen::VectorXd unitTangent(numberOfCoords_);

    for (int i = 0; i < numberOfStates_; ++i) {
        gradientCallCount_++;

        wf_.evaluate(x.segment(i*numberOfCoords_, numberOfCoords_));
        grad.segment(i*numberOfCoords_, numberOfCoords_) =
                wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();

        stateGrad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();

        unitTangent = unitTangents_.row(i);
        grad.segment(i*numberOfCoords_, numberOfCoords_) = stateGrad - (stateGrad.dot(unitTangent))*unitTangent;

    }

    //std::cout << "grad\n" << grad.transpose() << std::endl;
    //std::cout << "grad orthogonal\n" << orthogonalGrad.transpose() << std::endl;
}

bool StringOptimizationProblem::callback(cppoptlib::Criteria<double> &state, Eigen::VectorXd &x) {
    stepCounter_++;
    std::cout << "(" << std::setw(2) << state.iterations << ")"
              << " f(x) = "     << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
              << " xDelta = "   << std::setw(8) << state.xDelta
              << std::endl;
    return true;
}

