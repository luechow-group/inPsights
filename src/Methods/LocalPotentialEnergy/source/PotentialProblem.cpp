/* Copyright (C) 2018 Leonard Reuter.
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

#include "PotentialProblem.h"
#include <ParticlesVector.h>
#include <Metrics.h>

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

    ElectronsVector ec(PositionsVector{x});

    // initiating V_pot as nn potentials
    value = vPotNuclei_;

    // calculating en potentials
    for ( int i = 0; i < ec.numberOfEntities(); i++) {
        for (int j = 0; j < nuclei_.numberOfEntities(); j++) {
            value += -1.0 * nuclei_[j].charge() \
                / Metrics::distance(ec[i].position(), nuclei_[j].position());
        }
    }

    // calculating ee potentials
    for ( int i = 0; i < ec.numberOfEntities(); i++) {
        for (int j = i + 1; j < ec.numberOfEntities(); j++) {
            value += 1.0 \
                / Metrics::distance(ec[i].position(), ec[j].position());
        }
    }

    return value;
}

void PotentialProblem::gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad){
    int ei; // electron index
    int ci; // coordinate index (0,1,2 =: x,y,z)
    double sum;

    ElectronsVector ec(PositionsVector{x});

    for (int i = 0; i < x.size(); i++){
        ei = i / 3;
        ci = i % 3;
        sum = 0.0;

        // calculating en contributions
        for (int j = 0; j < nuclei_.numberOfEntities(); j++){
            sum += nuclei_[j].charge() \
                * (nuclei_[j].position()[ci] \
                    - ec[ei].position()[ci]) \
                / pow(Metrics::distance(ec[ei].position(), nuclei_[j].position()), 3.0);
        }

        // calculating ee contributions
        for (int j = 0; j < ec.numberOfEntities(); j++){
            if (j != ei){
                sum += -1.0 \
                * (ec[j].position()[ci] \
                    - ec[ei].position()[ci]) \
                / pow(Metrics::distance(ec[ei].position(), ec[j].position()), 3.0);
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
                / Metrics::distance(nuclei_[i].position(), nuclei_[j].position());
        }
    }
}

