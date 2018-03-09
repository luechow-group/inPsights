//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_SPINTYPECOLLECTION_H
#define AMOLQCGUI_SPINTYPECOLLECTION_H

#include <Eigen/Core>
#include "SpinType.h"
#include "AbstractCollection.h"

class SpinTypeCollection : public AbstractCollection{
public:
    explicit SpinTypeCollection(long size = 0);
    SpinTypeCollection(unsigned long numberOfAlphaElectrons, unsigned long numberOfBetaElectrons);

    explicit SpinTypeCollection(const Eigen::VectorXi& spinTypes);

    Spin::SpinType operator[](long i) const;

    void insert(Spin::SpinType spinType, long i);
    void append(Spin::SpinType spinType);
    void prepend(Spin::SpinType spinType);
    void permute(long i, long j) override;

    const Eigen::VectorXi& spinTypesAsEigenVector() const;

    Eigen::VectorXi& spinTypesAsEigenVector();

    friend std::ostream& operator<<(std::ostream& os, const SpinTypeCollection& pc);

private:
    Eigen::VectorXi spinTypes_;
};

#endif //AMOLQCGUI_SPINTYPECOLLECTION_H
