//
// Created by Michael Heuer on 13.01.17.
//

#ifndef AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
#define AMOLQCGUI_ELECTRONICWAVEFUNCTION_H

#include <Eigen/Core>
#include <vector>
#include "AtomCollection.h"
#include "ElectronCollections.h"

namespace ElectronPositioningMode {
  typedef enum {
    GAUSSIAN, DENSITY, LMO
  } electronPositioningModeType;
}

extern "C" {
void amolqc_init();
void amolqc_set_wf(int *nElecs, int *nAtoms, const char* fileName);
void amolqc_initial_positions(ElectronPositioningMode::electronPositioningModeType mode, int nElecs, double x[]);
void amolqc_eloc(double x[], int n, double *phi, double *u, double grad[], double *elocal);
}

class ElectronicWaveFunction {

public:
    static ElectronicWaveFunction& getInstance(const std::string& fileName = "");

    const std::string& getFileName();

    void initialize(const std::string& fileName);

    void setRandomElectronPositionCollection(unsigned electronNumber,
                                           ElectronPositioningMode::electronPositioningModeType);

    void evaluate(const ElectronCollection& electronCollection);

    void evaluate(const Eigen::VectorXd &electronPositionCollection);

    double getLocalEnergy();

    double getDeterminantProbabilityAmplitude();

    double getJastrowFactor();

    double getProbabilityAmplitude();

    double getProbabilityDensity();

    double getNegativeLogarithmizedProbabilityDensity();

    double getInverseNegativeLogarithmizedProbabilityDensity();

    Eigen::VectorXd getElectronPositionCollection();

    Eigen::VectorXd getElectronDriftCollection();

    Eigen::VectorXd getProbabilityAmplitudeGradientCollection();

    Eigen::VectorXd getProbabilityDensityGradientCollection();

    Eigen::VectorXd getNegativeLogarithmizedProbabilityDensityGradientCollection();

    Eigen::VectorXd getInverseNegativeLogarithmizedProbabilityDensityGradientCollection();

    unsigned long getNumberOfNuclei() const;

    AtomCollection getAtomCollection();

    unsigned long getNumberOfElectrons() const;

    SpinTypeCollection getSpinTypeCollection() const;

    void setFrozenElectrons(const std::vector<unsigned long>& frozenElectrons);
    std::vector<unsigned long> getFrozenElectrons();

private:
    explicit ElectronicWaveFunction(const std::string& fileName);
    const std::string fileName_;
    unsigned long numberOfNuclei_, numberOfElectrons_, numberOfAlphaElectrons_, numberOfBetaElectrons_;
    double determinantProbabilityAmplitude_, jastrowFactor_, localEnergy_;
    Eigen::VectorXd electronPositionCollection_, electronDriftCollection_;

    std::vector<unsigned long> frozenElectrons_;
    AtomCollection atomCollection_;
    SpinTypeCollection spinTypeCollection_;

};

#endif //AMOLQCGUI_ELECTRONICWAVEFUNCTION_H
