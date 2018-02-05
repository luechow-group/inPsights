//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_SPINTYPECOLLECTION_H
#define AMOLQCGUI_SPINTYPECOLLECTION_H

#include <Eigen/Core>
#include "SpinType.h"

class SpinTypeCollection{
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
    unsigned long numberOfSpinTypes_;
    Eigen::VectorXi spinTypes_;
};

#endif //AMOLQCGUI_SPINTYPECOLLECTION_H
