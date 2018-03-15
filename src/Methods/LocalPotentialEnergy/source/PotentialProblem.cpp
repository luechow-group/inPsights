//
// Created by leonard on 22.02.18.
//

#include "PotentialProblem.h"
#include "ElectronsVector.h"

PotentialProblem::PotentialProblem()
        : vPotNuclei_(0) {
    calculateVPotNuclei();
}

PotentialProblem::PotentialProblem(const AtomsVector & atomsVector)
        : nuclei_(atomsVector),
          vPotNuclei_(0) {
    calculateVPotNuclei();
}

double PotentialProblem::value(const Eigen::VectorXd &x) {
    double value;

    ElectronsVector ec(x);

    // initiating V_pot as nn potentials
    value = vPotNuclei_;

    // calculating en potentials
    for ( int i = 0; i < ec.numberOfEntities(); i++) {
        for (int j = 0; j < nuclei_.numberOfEntities(); j++) {
            value += -1.0 * nuclei_[j].charge() \
                / Particle::distance(ec[i], nuclei_[j]);
        }
    }

    // calculating ee potentials
    for ( int i = 0; i < ec.numberOfEntities(); i++) {
        for (int j = i + 1; j < ec.numberOfEntities(); j++) {
            value += 1.0 \
                / Particle::distance(ec[i], ec[j]);
        }
    }

    return value;
}

void PotentialProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad){
    int ei; // electron index
    int ci; // coordinate index (0,1,2 =: x,y,z)
    double sum;

    ElectronsVector ec(x);

    for (int i = 0; i < x.size(); i++){
        ei = i / 3;
        ci = i % 3;
        sum = 0.0;

        // calculating en contributions
        for (int j = 0; j < nuclei_.numberOfEntities(); j++){
            sum += nuclei_[j].charge() \
                * (nuclei_[j].position()[ci] \
                    - ec[ei].position()[ci]) \
                / pow(Particle::distance(ec[ei], nuclei_[j]), 3.0);
        }

        // calculating ee contributions
        for (int j = 0; j < ec.numberOfEntities(); j++){
            if (j != ei){
                sum += -1.0 \
                * (ec[j].position()[ci] \
                    - ec[ei].position()[ci]) \
                / pow(Particle::distance(ec[ei], ec[j]), 3.0);
            }
        }

        grad(i) = -1.0 * sum;
    }
}

void PotentialProblem::hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &hessian) {
}


AtomsVector PotentialProblem::getAtoms() const{
    return nuclei_;
}

void PotentialProblem::calculateVPotNuclei() {
    vPotNuclei_ = 0.0;
    for (int i = 0; i < nuclei_.numberOfEntities(); i++){
        for (int j = i+1; j< nuclei_.numberOfEntities(); j++){
            vPotNuclei_ += nuclei_[i].charge() * nuclei_[j].charge() \
                / Particle::distance(nuclei_[i], nuclei_[j]);
        }
    }
}

