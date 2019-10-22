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

#ifndef INPSIGHTS_ELECTRONICWAVEFUNCTION_H
#define INPSIGHTS_ELECTRONICWAVEFUNCTION_H

#include <Eigen/Core>
#include <vector>

#include <ParticlesVectorCollection.h>

namespace ElectronPositioningMode {
    typedef enum {
        GAUSSIAN, DENSITY, LMO
    } electronPositioningModeType;
}

extern "C" {
void amolqc_init();
void amolqc_set_wf(int *nElecs, int *nAtoms, const char *fileName);
void amolqc_initial_positions(ElectronPositioningMode::electronPositioningModeType mode, int nElecs, double x[]);
void amolqc_eloc(double x[], int n, double *phi, double *u, double grad[], double *elocal);
}

class ElectronicWaveFunction {

public:
    static ElectronicWaveFunction &getEmpty();

    static ElectronicWaveFunction &getInstance(const std::string &fileName = "");

    const std::string &getFileName();

    void initialize(const std::string &fileName);

    void
    setRandomElectronPositionsVector(unsigned electronNumber, ElectronPositioningMode::electronPositioningModeType);

    void evaluate(const ElectronsVector &electronsVector);

    void evaluate(const Eigen::VectorXd &electronPositionsVector);

    double getLocalEnergy();

    double getDeterminantProbabilityAmplitude();

    double getJastrowFactor();

    double getProbabilityAmplitude();

    double getProbabilityDensity();

    double getNegativeLogarithmizedProbabilityDensity();

    double getInverseNegativeLogarithmizedProbabilityDensity();

    ElectronsVector getElectronsVector();

    Eigen::VectorXd getElectronDriftCollection();

    Eigen::VectorXd getProbabilityAmplitudeGradientCollection();

    Eigen::VectorXd getProbabilityDensityGradientCollection();

    Eigen::VectorXd getNegativeLogarithmizedProbabilityDensityGradientCollection();

    Eigen::VectorXd getInverseNegativeLogarithmizedProbabilityDensityGradientCollection();

    unsigned long getNumberOfNuclei() const;

    AtomsVector getAtomsVector() const;

    unsigned long getNumberOfElectrons() const;

    SpinTypesVector getSpinTypesVector() const;

private:
    explicit ElectronicWaveFunction();

    explicit ElectronicWaveFunction(const std::string &fileName);

    const std::string fileName_;
    unsigned long numberOfNuclei_, numberOfElectrons_, numberOfAlphaElectrons_, numberOfBetaElectrons_;
    double determinantProbabilityAmplitude_, jastrowFactor_, localEnergy_;
    Eigen::VectorXd electronPositionsVectorAsEigenVector_, electronDriftCollection_;//TODO REPLACE BY BASELIB ELECTRONCOLLECTION!
    AtomsVector atomsVector_;
    SpinTypesVector spinTypesVector_;

};

#endif //INPSIGHTS_ELECTRONICWAVEFUNCTION_H
