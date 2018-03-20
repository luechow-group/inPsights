//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_ELECTRONCOLLECTION_H
#define AMOLQCPP_ELECTRONCOLLECTION_H

#include "ParticlesVector.h"
#include "SpinTypesVector.h"
#include "Electron.h"

class ElectronsVector : public ParticlesVector{
public:
    ElectronsVector() = default; //TODO check!!
    explicit ElectronsVector(const Eigen::VectorXd& positions);
    ElectronsVector(const Eigen::VectorXd& positions, const Eigen::VectorXi& spinTypes);

    ElectronsVector(const PositionsVector& positionsVector,
                       const SpinTypesVector& spinTypesVector);

    Electron operator[](long i) const;

    void insert (const Electron& electron, long i);
    void append (const Electron& electron);
    void prepend(const Electron& electron);
    void permute(long i, long j) override;

    friend std::ostream& operator<<(std::ostream& os, const ElectronsVector& ec);

    const SpinTypesVector& spinTypesVector() const;
    SpinTypesVector& spinTypesVector();

private:
    SpinTypesVector spinTypesVector_;
};

#endif //AMOLQCPP_ELECTRONCOLLECTION_H
