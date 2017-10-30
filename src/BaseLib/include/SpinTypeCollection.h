//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_SPINTYPECOLLECTION_H
#define AMOLQCGUI_SPINTYPECOLLECTION_H

#include <Eigen/Core>
#include "SpinType.h"

class SpinTypeCollection{
public:
    explicit SpinTypeCollection(long size = 0);
    explicit SpinTypeCollection(const Eigen::VectorXi& spinTypes);

    Spin::SpinType spinType(long i);

    unsigned long numberOfSpinTypes() const;

    void insert(Spin::SpinType spinType, long i);
    void append(Spin::SpinType spinType);
    void prepend(Spin::SpinType spinType);

    void setSpinType(long i, Spin::SpinType spinType);

    Eigen::VectorXi spinTypesAsEigenVector();

private:
    unsigned long numberOfSpinTypes_;
    Eigen::VectorXi spinTypes_;
};

#endif //AMOLQCGUI_SPINTYPECOLLECTION_H
