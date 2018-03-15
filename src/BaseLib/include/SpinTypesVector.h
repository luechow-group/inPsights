//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_SPINTYPECOLLECTION_H
#define AMOLQCPP_SPINTYPECOLLECTION_H

#include <Eigen/Core>
#include "SpinType.h"
#include "AbstractVector.h"

class SpinTypesVector : public AbstractVector{
public:
    explicit SpinTypesVector(long size = 0);
    SpinTypesVector(unsigned long numberOfAlphaElectrons, unsigned long numberOfBetaElectrons);

    explicit SpinTypesVector(const Eigen::VectorXi& spinTypes);

    Spin::SpinType operator[](long i) const;

    void insert(Spin::SpinType spinType, long i);
    void append(Spin::SpinType spinType);
    void prepend(Spin::SpinType spinType);
    void permute(long i, long j) override;

    const Eigen::VectorXi& spinTypesAsEigenVector() const;

    Eigen::VectorXi& spinTypesAsEigenVector();

    friend std::ostream& operator<<(std::ostream& os, const SpinTypesVector& pc);

private:
    Eigen::VectorXi spinTypes_;
};

#endif //AMOLQCPP_SPINTYPECOLLECTION_H
