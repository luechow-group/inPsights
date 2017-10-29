//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_SPINTYPECOLLECTION_H
#define AMOLQCGUI_SPINTYPECOLLECTION_H

#include <Eigen/Core>
#include "SpinType.h"

using namespace Eigen;



class SpinTypeCollection{
public:
    explicit SpinTypeCollection(long size);

    explicit SpinTypeCollection(const VectorXi& spinTypes);

    Spin::SpinType spinType(long i);

    void setSpinType(long i, Spin::SpinType spinType);

    VectorXi asVectorXi();

private:
    long size_;
    VectorXi spinTypes_;
};

#endif //AMOLQCGUI_SPINTYPECOLLECTION_H
