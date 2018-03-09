//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTION_H
#define AMOLQCGUI_ELECTRONCOLLECTION_H

#include "ParticleCollection.h"
#include "SpinTypeCollection.h"
#include "Electron.h"

class ElectronCollection : public ParticleCollection{
public:
    ElectronCollection() = default; //TODO check!!
    explicit ElectronCollection(const Eigen::VectorXd& positions);
    ElectronCollection(const Eigen::VectorXd& positions, const Eigen::VectorXi& spinTypes);

    ElectronCollection(const PositionsVector& positionsVector,
                       const SpinTypeCollection& spinTypeCollection);

    Electron operator[](long i) const;

    void insert (const Electron& electron, long i);
    void append (const Electron& electron);
    void prepend(const Electron& electron);
    void permute(long i, long j) override;

    friend std::ostream& operator<<(std::ostream& os, const ElectronCollection& ec);

    const SpinTypeCollection& spinTypeCollection() const;
    SpinTypeCollection& spinTypeCollection();

private:
    SpinTypeCollection spinTypeCollection_;
};

#endif //AMOLQCGUI_ELECTRONCOLLECTION_H
