//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTION_H
#define AMOLQCGUI_ELECTRONCOLLECTION_H

#include "ParticleCollection.h"

#include "Electron.h"

class ElectronCollection : public ParticleCollection{
public:
    ElectronCollection() = default;
    explicit ElectronCollection(const Eigen::VectorXd& positions);
    ElectronCollection(const Eigen::VectorXd& positions, const Eigen::VectorXi& spinTypes);
    ElectronCollection(const ParticleCollection& particleCollection, const Eigen::VectorXi& spinTypes);

    Electron electron(long i) const;

    void insert (const Electron& electron, long i);
    void append (const Electron& electron);
    void prepend(const Electron& electron);
    void permute(long i, long j) override ;

    friend std::ostream& operator<<(std::ostream& os, const ElectronCollection& ec);

    Eigen::VectorXi getSpinTypes() const;

private:
    Eigen::VectorXi spinTypes_;
};

#endif //AMOLQCGUI_ELECTRONCOLLECTION_H


/*
#include <Eigen/Core>
#include "SpinType.h"

class SpinTypeCollection : public AbstractCollection{
public:
    explicit SpinTypeCollection(unsigned long size = 0);
    SpinTypeCollection(unsigned long numberOfAlphaElectrons, unsigned long numberOfBetaElectrons);

    explicit SpinTypeCollection(const Eigen::VectorXi& spinTypes);

    Spin::SpinType spinType(long i) const;

    unsigned long numberOfSpinTypes() const;

    void insert(Spin::SpinType spinType, long i);
    void append(Spin::SpinType spinType);
    void prepend(Spin::SpinType spinType);
    void permute(long i, long j);

    void setSpinType(long i, Spin::SpinType spinType);

    Eigen::VectorXi spinTypesAsEigenVector() const;

    friend std::ostream& operator<<(std::ostream& os, const SpinTypeCollection& pc);

private:
    Eigen::VectorXi spinTypes_;
};
 */