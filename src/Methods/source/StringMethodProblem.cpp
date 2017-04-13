//
// Created by heuer on 13.04.17.
//

#include <meta.h>
#include "StringMethodProblem.h"
#include "PointInterpolationGenerator.h"

StringMethodProblem::StringMethodProblem(long numberOfStates,
                                         long numberOfCoords)
        : numberOfStates_(numberOfStates),
          numberOfCoords_(numberOfCoords),
          valueCallCount_(0),
          gradientCallCount_(0)
{
    assert(numberOfStates_ > 2);
    assert( (numberOfCoords_ > 0) && ( numberOfCoords_%3 == 0) );
    chain_.resize(numberOfStates,numberOfCoords+1);
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
    for (int i = 0; i < numberOfStates_; ++i) {
        valueCallCount_++;
        Eigen::VectorXd stateCoordVector(systemCoordVector.segment(i*numberOfCoords_, numberOfCoords_) );

        wf_.evaluate(stateCoordVector);
        chain_(i,0) = wf_.getNegativeLogarithmizedProbabilityDensity();
        //chain_.block(i,1,1,numberOfCoords_) = stateCoordVector;
    }

    return chain_.col(0).sum();
}

void StringMethodProblem::gradient(const Eigen::VectorXd &systemCoordVector,
                                   Eigen::VectorXd &systenGradientVector) {


    for (int i = 0; i < numberOfStates_; ++i) {
        gradientCallCount_++;

        wf_.evaluate(systemCoordVector.segment(i*numberOfCoords_, numberOfCoords_));
        chainGradient_.segment(i*numberOfCoords_, numberOfCoords_) =
                wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
    }

    // calculate the single state gradients orthogonal to string

    // make one giant system state vector out of
}

bool StringMethodProblem::callback(const cppoptlib::Criteria<double> &state, const Eigen::VectorXd &x) {
    return false;
}

void StringMethodProblem::resetString() {

}

void StringMethodProblem::reparametrizeString() {
    // include probability density in the interpolation

    // construct (value weighted) arclength parametrized spline

    BSplines::PointInterpolationGenerator interpolator(chain_,3);
    /*Eigen::VectorXi excludedDimensions(1);
    excludedDimensions << 1;
    arcLengthParametrizedBSpline_ = BSplines::ArcLengthParametrizedBSpline(interpolator.generateBSpline(1),
                                                                           excludedDimensions);*/
    bSpline_ = interpolator.generateBSpline(1);
}


void StringMethodProblem::constructChainOfStates() {



    // how many states?
    //wf_.evaluate(x);
}
