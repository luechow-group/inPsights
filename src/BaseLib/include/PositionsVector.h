//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_POSITIONCOLLECTION_H
#define AMOLQCPP_POSITIONCOLLECTION_H

#include <Eigen/Core>
#include "AbstractVector.h"

class PositionsVector : public AbstractVector{
public:
    PositionsVector();
    explicit PositionsVector(const Eigen::VectorXd& positions);

    Eigen::Vector3d operator[](long i) const;

    void insert(const Eigen::Vector3d& position, long i);
    void append(const Eigen::Vector3d& position);
    void prepend(const Eigen::Vector3d& position);
    void permute(long i, long j) override;

    /* TODO
    void replace(long i);
    void remove(long i);
    ParticleCollection part(std::vector<long> indices);
    */

    friend std::ostream& operator<<(std::ostream& os, const PositionsVector& pc);

    const Eigen::VectorXd & positionsAsEigenVector() const;

    Eigen::VectorXd & positionsAsEigenVector();

private:
    Eigen::VectorXd positions_;
    unsigned entityLength_ = 3;

    long calculateIndex(long i) const override ;
};

#endif //AMOLQCPP_POSITIONCOLLECTION_H
