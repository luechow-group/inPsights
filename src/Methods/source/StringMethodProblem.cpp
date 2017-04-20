//
// Created by heuer on 13.04.17.
//

#include <meta.h>
#include "StringMethodProblem.h"
#include "PointInterpolationGenerator.h"

StringMethodProblem::StringMethodProblem(long numberOfStates,
                                         long numberOfCoords)
        : stepCounter_(0),
          numberOfStates_(numberOfStates),
          numberOfCoords_(numberOfCoords),
          chain_(numberOfStates,numberOfCoords),
          //chainGradient_(numberOfStates*numberOfCoords),
          unitTangentVector_(numberOfStates*numberOfCoords),
          uValues_(numberOfStates),
          valueCallCount_(0),
          gradientCallCount_(0)
{
    assert(numberOfStates_ > 2);
    assert( (numberOfCoords_ > 0) && ( numberOfCoords_%3 == 0) );
}
/*StringMethodProblem::StringMethodProblem(const Eigen::MatrixXd &coords)
        : chain_(Eigen::MatrixXd(coords.rows(),coords.cols()+1)),
          chainGradient_(Eigen::VectorXd(coords.rows()*coords.cols())),
          numberOfStates_(chain_.rows()),
          numberOfCoords_(chain_.cols()-1),
          valueCallCount_(0),
          gradientCallCount_(0)
{
    assert(numberOfStates_ > 2);
    assert( (numberOfCoords_ > 0) && ( numberOfCoords_%3 == 0) );

    chain_.col(0) = Eigen::VectorXd::Zero(numberOfStates_);
    chain_.rightCols(numberOfCoords_) = coords;
}*/


double StringMethodProblem::value(const Eigen::VectorXd &systemCoordVector) {

    double value = 0;
    for (int i = 0; i < numberOfStates_; ++i) {
        Eigen::VectorXd stateCoordVector(systemCoordVector.segment(i*numberOfCoords_, numberOfCoords_) );

        valueCallCount_++;
        wf_.evaluate(stateCoordVector);

        value += wf_.getNegativeLogarithmizedProbabilityDensity();
        //chain_(i,0) = wf_.getNegativeLogarithmizedProbabilityDensity();
    }
    return value;//chain_.col(0).sum();
}

void StringMethodProblem::gradient(const Eigen::VectorXd &systemCoordVector,
                                   Eigen::VectorXd &systemGradientVector) {
    Eigen::VectorXd chainGradient_(systemCoordVector.size());
    for (int i = 0; i < numberOfStates_; ++i) {
        gradientCallCount_++;

        wf_.evaluate(systemCoordVector.segment(i*numberOfCoords_, numberOfCoords_));
        chainGradient_.segment(i*numberOfCoords_, numberOfCoords_) =
                wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
    }

    calculateUnitTangentVector();
    systemGradientVector = chainGradient_ - (chainGradient_.dot(unitTangentVector_))*unitTangentVector_;
}

bool StringMethodProblem::callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x) {
    stepCounter_++;
    std::cout << "callback" << std::endl;
    setChainOfStates(x);

    if(stepCounter_ > 4) {
        stepCounter_ = 0;
        reparametrizeString();
        std::cout << "reparametrized" << std::endl;
    }

    return true;
}

void StringMethodProblem::getChainOfStatesFromSpline() {
    double u;

    for (int i = 0; i < numberOfStates_ ; ++i) {
        u = double(i)/double(numberOfStates_-1);
        uValues_(i) = u;
        chain_.row(i) = bSpline_.evaluate(u);
        //systemGradientVector.segment(i*numberOfCoords_, numberOfCoords_) =
        //unitTangentVector_.block(i*numberOfCoords_,1,numberOfCoords_,1) = bSpline_.evaluate(uValues_(i),1);
    }
}

void StringMethodProblem::setChainOfStates(const Eigen::VectorXd &systemCoordVector) {
    for (int i = 0; i < numberOfStates_; ++i) {
        chain_.row(i) = systemCoordVector.segment(i*numberOfCoords_, numberOfCoords_);
    }
}

void StringMethodProblem::calculateUnitTangentVector() {

    for (int i = 0; i < numberOfStates_ ; ++i) {
        unitTangentVector_.block(i*numberOfCoords_,0,numberOfCoords_,1) = bSpline_.evaluate(uValues_(i),1);
    }
    unitTangentVector_.normalize();
}


void StringMethodProblem::resetOptimizer() {

}

void StringMethodProblem::reparametrizeString() {
    // include probability density in the interpolation
    // construct (value weighted) arclength parametrized spline

    BSplines::PointInterpolationGenerator interpolator(chain_,3,true);
    /*Eigen::VectorXi excludedDimensions(1);
    excludedDimensions << 1;
    arcLengthParametrizedBSpline_ = BSplines::ArcLengthParametrizedBSpline(interpolator.generateBSpline(1),
                                                                           excludedDimensions);*/
    bSpline_ = interpolator.generateBSpline(1);
}
