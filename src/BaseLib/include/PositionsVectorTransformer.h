//
// Created by Michael Heuer on 08.04.18.
//

#ifndef AMOLQCPP_POSITIONSVECTORTRANSFORMER_H
#define AMOLQCPP_POSITIONSVECTORTRANSFORMER_H

#include "PositionsVector.h"
#include <Eigen/Eigenvalues>
#include <vector>

namespace PositionsVectorTransformer{

    struct AngleAxis{
        AngleAxis(double angle,Eigen::Vector3d axis)
                : angle_(angle), axis_(axis){};
        double angle_;
        Eigen::Vector3d axis_;
    };

    void translateCenterOfMassToOrigin(PositionsVector& positionsVector);

    void rotateAroundAxis(PositionsVector &p, double angle,
                          const Eigen::Vector3d &axisStart,
                          const Eigen::Vector3d &axisEnd);
    void rotateAroundAxis(PositionsVector &positionsVector, double angle,
                          const Eigen::Vector3d &axis);

    Eigen::Matrix3d rotationMatrixFromQuaternion(const Eigen::Vector4d &q);

    Eigen::Vector3d calculateCenterOfMass(const PositionsVector& positionsVector,
                                   const Eigen::VectorXd& weights);

    Eigen::Vector3d calculateCenterOfMass(const PositionsVector& positionsVector);

    Eigen::Vector4d quaternionFromAngleAndAxis(double angle, const Eigen::Vector3d &axis);

    AngleAxis quaternionToAngleAndAxis(const Eigen::Vector4d &quaternion);
};


#endif //AMOLQCPP_POSITIONSVECTORTRANSFORMER_H