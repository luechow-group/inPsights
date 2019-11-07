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

#include <ElectronicWaveFunctionProblem.h>
#include <iomanip>

ElectronicWaveFunctionProblem::ElectronicWaveFunctionProblem()
        :
        valueCallCount_(0),
        gradientCallCount_(0),
        wf_(ElectronicWaveFunction::getEmpty()),
        optimizationPath_(wf_.getSpinTypesVector()),
        electronCoordinateIndicesThatWereNaN_(Eigen::VectorXb(wf_.getNumberOfElectrons()*3).setConstant(false)),
        indicesOfElectronsNotAtNuclei_(0),
        indicesOfElectronsAtNuclei_(0)
{}

ElectronicWaveFunctionProblem::ElectronicWaveFunctionProblem(const std::string &fileName, const bool &putElectronsIntoNuclei, const bool &printStatus)
        :
        valueCallCount_(0),
        gradientCallCount_(0),
        wf_(ElectronicWaveFunction::getInstance(fileName)),
        optimizationPath_(wf_.getSpinTypesVector()),
        electronCoordinateIndicesThatWereNaN_(Eigen::VectorXb(wf_.getNumberOfElectrons()*3).setConstant(false)),
        indicesOfElectronsNotAtNuclei_(0),
        indicesOfElectronsAtNuclei_(0),
        putElectronsIntoNuclei_(putElectronsIntoNuclei),
        printStatus_(printStatus)
{
    for (unsigned long i = 0; i < wf_.getNumberOfElectrons(); ++i) {
        indicesOfElectronsNotAtNuclei_.push_back(i);
    }
}

double ElectronicWaveFunctionProblem::value(const Eigen::VectorXd &x) {
    valueCallCount_++;
    wf_.evaluate(x);

    return wf_.getNegativeLogarithmizedProbabilityDensity();
}

void ElectronicWaveFunctionProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    gradientCallCount_++;
    wf_.evaluate(x);

    grad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
    fixGradient(grad);
}

void ElectronicWaveFunctionProblem::hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &hessian) {

    //TODO asserts?
    long dims = x.size();

    cppoptlib::Problem<double,Eigen::Dynamic>::semifiniteHessian(x,hessian,indicesOfElectronsNotAtNuclei_,3,0);
}

void ElectronicWaveFunctionProblem::fixGradient(Eigen::VectorXd &gradient) {
    for(auto i : indicesOfElectronsNotAtNuclei_){
        for (unsigned j = 0; j < 3; ++j) {
            if(gradient[i*3+j] != gradient[i*3+j]) {
                gradient[i*3+j] = 0;
                electronCoordinateIndicesThatWereNaN_[i*3+j] = true;
            }
            else electronCoordinateIndicesThatWereNaN_[i*3+j] = false;
        }
    }
    for(auto i : indicesOfElectronsAtNuclei_){
        gradient.segment(i*3,3) = Eigen::Vector3d::Zero(3);
    }
}


void ElectronicWaveFunctionProblem::putElectronsIntoNuclei(Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    assert( x.size() == wf_.getNumberOfElectrons()*3 && "Number of dimensions must be identical and multiple of 3");

    auto atomsVector = wf_.getAtomsVector();
    auto numberOfNuclei = atomsVector.numberOfEntities();
    auto numberOfElectrons = wf_.getNumberOfElectrons();

    // iterate over electrons that were not at nuclei in the last step
    for(auto i : indicesOfElectronsNotAtNuclei_){
        unsigned long closestNucleusIdx=0;
        double smallestDistance = std::numeric_limits<double>::infinity();

        // iterate over all nuclei and find the index of the electron with the smallest distance
        for(unsigned long j = 0; j < numberOfNuclei; ++j){
            double distance = (atomsVector[j].position()-x.segment(i*3,3)).norm();
            if(distance < smallestDistance ) {
                smallestDistance = distance;
                closestNucleusIdx = j;
            }
        }
        // check the electron with the smallest distance is close than the threshold
        double threshold = 0.0005;
        if (smallestDistance <= threshold){
            //TODO PROPER?
            //Eigen::Block<Eigen::VectorXd, i*3, 0>(x.derived(), 0, 0) = atomsVector[closestNucleusIdx].position();
            x.segment(i*3,3) = atomsVector[closestNucleusIdx].position();

            //save the electron index in the indicesOfElectronsAtNuclei_ vector
            indicesOfElectronsAtNuclei_.push_back(i);
            std::sort(indicesOfElectronsAtNuclei_.begin(),indicesOfElectronsAtNuclei_.end());

            // recalulate gradient
            gradient(x,grad);
            gradientResetQ = true;
        }
    }
    // now, remove the electrons from the indicesOfElectronsNotAtNuclei_ vector
    for(auto i : indicesOfElectronsAtNuclei_) {
        indicesOfElectronsNotAtNuclei_.erase(std::remove(indicesOfElectronsNotAtNuclei_.begin(),
                                                         indicesOfElectronsNotAtNuclei_.end(), i),
                                             indicesOfElectronsNotAtNuclei_.end());
    }
}

AtomsVector ElectronicWaveFunctionProblem::getAtomsVector() const{
    return wf_.getAtomsVector();
}

std::vector<unsigned long> ElectronicWaveFunctionProblem::getIndicesOfElectronsNotAtNuclei() {
    return indicesOfElectronsNotAtNuclei_;
}

std::vector<unsigned long> ElectronicWaveFunctionProblem::getIndicesOfElectronsAtNuclei() {
    return indicesOfElectronsAtNuclei_;
}