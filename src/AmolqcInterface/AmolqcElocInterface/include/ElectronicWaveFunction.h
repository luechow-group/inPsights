// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
