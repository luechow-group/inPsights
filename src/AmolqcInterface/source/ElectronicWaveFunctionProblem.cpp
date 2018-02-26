//
// Created by heuer on 06.04.17.
//

#include "ElectronicWaveFunctionProblem.h"
#include <iomanip>

ElectronicWaveFunctionProblem::ElectronicWaveFunctionProblem(const std::string &fileName)
        :
        valueCallCount_(0),
        gradientCallCount_(0),
        wf_(ElectronicWaveFunction::getInstance(fileName)),
        optimizationPath_(wf_.getSpinTypeCollection()),
        electronCoordinateIndicesThatWereNaN_(Eigen::Matrix<bool,Eigen::Dynamic,1>(wf_.getNumberOfElectrons()*3).setConstant(false)),
        indicesOfElectronsNotAtNuclei_(0),
        indicesOfElectronsAtNuclei_(0)
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

    cppoptlib::Problem<double,Eigen::Dynamic>::semifiniteHessian(x,hessian,2);

    for (auto i : indicesOfElectronsAtNuclei_){
        hessian.block(i*3,0,3,dims) = Eigen::MatrixXd::Zero(3,dims);
        hessian.block(0,i*3,dims,3) = Eigen::MatrixXd::Zero(dims,3);
    }
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

    auto atomCollection = wf_.getAtomCollection();
    auto numberOfNuclei = atomCollection.numberOfParticles();
    auto numberOfElectrons = wf_.getNumberOfElectrons();

    // iterate over electrons that were not at nuclei in the last step
    for(auto i : indicesOfElectronsNotAtNuclei_){
        unsigned long closestNucleusIdx=0;
        double smallestDistance = std::numeric_limits<double>::infinity();

        // iterate over all nuclei and find the index of the electron with the smallest distance
        for(unsigned long j = 0; j < numberOfNuclei; ++j){
            double distance = (atomCollection[j].position()-x.segment(i*3,3)).norm();
            if(distance < smallestDistance ) {
                smallestDistance = distance;
                closestNucleusIdx = j;
            }
        }
        // check the electron with the smallest distance is close than the threshold
        double threshold = 0.0005;
        if (smallestDistance <= threshold){
            //TODO PROPER?
            //Eigen::Block<Eigen::VectorXd, i*3, 0>(x.derived(), 0, 0) = atomCollection[closestNucleusIdx].position();
            x.segment(i*3,3) = atomCollection[closestNucleusIdx].position();

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

bool ElectronicWaveFunctionProblem::callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) {
    gradientResetQ = false;
    putElectronsIntoNuclei(x, grad); //gradientQ could be true now

    optimizationPath_.append(ElectronCollection(x, wf_.getSpinTypeCollection().spinTypesAsEigenVector()));

    std::cout << "(" << std::setw(2) << state.iterations << ")"
              << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
              << " xDelta = " << std::setw(8) << state.xDelta
              << " gradInfNorm = " << std::setw(8) << state.gradNorm
              << std::endl;
    std::cout << "value calls: " <<  valueCallCount_ << ", gradient calls:" << gradientCallCount_ << std::endl;

    for (auto & it : indicesOfElectronsNotAtNuclei_) std::cout << it << " ";
    std::cout << std::endl;
    for (auto & it : indicesOfElectronsAtNuclei_) std::cout << it << " ";
    std::cout << std::endl;

    return true;
}

AtomCollection ElectronicWaveFunctionProblem::getAtomCollection() const{
    return wf_.getAtomCollection();
}

std::vector<unsigned long> ElectronicWaveFunctionProblem::getIndicesOfElectronsNotAtNuclei() {
    return indicesOfElectronsNotAtNuclei_;
}

std::vector<unsigned long> ElectronicWaveFunctionProblem::getIndicesOfElectronsAtNuclei() {
    return indicesOfElectronsAtNuclei_;
}
