// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_POSITIONSVECTORTRANSFORMER_H
#define INPSIGHTS_POSITIONSVECTORTRANSFORMER_H

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

    Eigen::Matrix3d rotationMatrixFromQuaternion(const Eigen::Vector4d &q);

    Eigen::Vector3d calculateCenterOfMass(const PositionsVector& positionsVector,
                                   const Eigen::VectorXd& weights);

    Eigen::Vector3d calculateCenterOfMass(const PositionsVector& positionsVector);

    Eigen::Vector4d quaternionFromAngleAndAxis(double angle, const Eigen::Vector3d &axis);

    AngleAxis quaternionToAngleAndAxis(const Eigen::Vector4d &quaternion);
};


#endif //INPSIGHTS_POSITIONSVECTORTRANSFORMER_H
