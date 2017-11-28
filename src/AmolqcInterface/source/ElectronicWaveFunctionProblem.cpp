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
        nanElectronsCoordinateIndices_(Eigen::Matrix<bool,Eigen::Dynamic,1>(wf_.getNumberOfElectrons()*3).setConstant(false)),
        indicesOfElectronsNotAtCoresYet_(0),
        indicesOfElectronsAtCores_(0)
{
    for (unsigned long i = 0; i < wf_.getNumberOfElectrons(); ++i) {
        indicesOfElectronsNotAtCoresYet_.push_back(i);
    }
}

double ElectronicWaveFunctionProblem::value(const Eigen::VectorXd &x) {
    valueCallCount_++;
    wf_.evaluate(x);
  return wf_.getNegativeLogarithmizedProbabilityDensity();
  //return wf_.getProbabilityDensity();
}

void ElectronicWaveFunctionProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
  gradientCallCount_++;
  wf_.evaluate(x);

    grad = wf_.getNegativeLogarithmizedProbabilityDensityGradientCollection();
    //std::cout << "unfixed " << grad.transpose() << std::endl;
    fixGradient(grad);
    //std::cout << "  fixed " << grad.transpose() << std::endl;
  //grad = wf_.getProbabilityDensityGradientCollection();
}

void ElectronicWaveFunctionProblem::hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &hessian) {

    //TODO asserts?
    long dims = x.size();
  //std::vector<unsigned long> ignoredIdx = wf_.getFrozenElectrons();
  //assert(ignoredIdx.size() < dims);

    cppoptlib::Problem<double,Eigen::Dynamic>::hessian(x,hessian);

    //TODO RETHINK!
    for (int i = 0; i < wf_.getNumberOfElectrons()*3; ++i) {
        if(nanElectronsCoordinateIndices_[i]){
            //set entire row to zero
            hessian.block(i,0,1,dims) = Eigen::MatrixXd::Zero(dims,1);
            //set entire column to zero
            hessian.block(0,i,dims,1) = Eigen::MatrixXd::Zero(1,dims);
        }
    }
}


void ElectronicWaveFunctionProblem::fixGradient(VectorXd &gradient) {

  //std::cout << "not at core: ";
    for(auto i : indicesOfElectronsNotAtCoresYet_){
  //    std::cout << i << " ";
        if(gradient[i*3+0] != gradient[i*3+0]) gradient[i*3+0] = 0;
        if(gradient[i*3+1] != gradient[i*3+1]) gradient[i*3+1] = 0;
        if(gradient[i*3+2] != gradient[i*3+2]) gradient[i*3+2] = 0;
    }
  //std::cout << std::endl;
  //std::cout << "    at core: ";
    for(auto i : indicesOfElectronsAtCores_){
  //    std::cout << i << " ";
        gradient.segment(i*3,3) = Eigen::Vector3d::Zero(3);
    }
  //std::cout << std::endl;
}


void ElectronicWaveFunctionProblem::putElectronsIntoNuclei(Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    assert( x.size() == wf_.getNumberOfElectrons()*3 && "Number of dimensions must be identical and multiple of 3");

    auto atomCollection = wf_.getAtomCollection();
    auto numberOfNuclei = atomCollection.numberOfParticles();
    auto numberOfElectrons = wf_.getNumberOfElectrons();

    for(auto i : indicesOfElectronsNotAtCoresYet_){
        unsigned long closestNucleusIdx=0;
        double smallestDistance = std::numeric_limits<double>::infinity();

        for(unsigned long j = 0; j < numberOfNuclei; ++j){
            double distance = (atomCollection[j].position()-x.segment(i*3,3)).norm();
            if(distance < smallestDistance ) {
                smallestDistance = distance;
                closestNucleusIdx = j;
            }
        }

        double threshold = 0.05;
        if (smallestDistance <= threshold){
            //TODO PROPER?
            //Eigen::Block<Eigen::VectorXd, i*3, 0>(x.derived(), 0, 0) = atomCollection[closestNucleusIdx].position();
            x.segment(i*3,3) = atomCollection[closestNucleusIdx].position();

            indicesOfElectronsNotAtCoresYet_.erase( std::remove( indicesOfElectronsNotAtCoresYet_.begin(),
                                                                 indicesOfElectronsNotAtCoresYet_.end(), i),
                                                    indicesOfElectronsNotAtCoresYet_.end() );
            indicesOfElectronsAtCores_.push_back(i);
            std::sort(indicesOfElectronsAtCores_.begin(),indicesOfElectronsAtCores_.end());

            // recalulate gradient
            gradient(x,grad);
        }

    }
}

bool ElectronicWaveFunctionProblem::callback(const cppoptlib::Criteria<double> &state, Eigen::VectorXd &x, Eigen::VectorXd& grad) {

    putElectronsIntoNuclei(x, grad);

    optimizationPath_.append(ElectronCollection(x, wf_.getSpinTypeCollection().spinTypesAsEigenVector()));

    std::cout << "(" << std::setw(2) << state.iterations << ")"
              << " f(x) = " << std::fixed << std::setw(8) << std::setprecision(8) << value(x)
              << " xDelta = " << std::setw(8) << state.xDelta
              << " gradInfNorm = " << std::setw(8) << state.gradNorm
              << std::endl;
    std::cout << "value calls: " <<  valueCallCount_ << ", gradient calls:" << gradientCallCount_ << std::endl;
    return true;
}

Eigen::VectorXd ElectronicWaveFunctionProblem::getNucleiPositions() const{
    return wf_.getAtomCollection().positionsAsEigenVector();
}

