//
// Created by Michael Heuer on 08.04.18.
//

#ifndef AMOLQCPP_POSITIONSVECTORTRANSFORMER_H
#define AMOLQCPP_POSITIONSVECTORTRANSFORMER_H

#include "PositionsVector.h"
#include <Eigen/Eigenvalues>
#include <vector>


//TODO refactor
Eigen::VectorXd permutePositionsCyclic(const Eigen::VectorXd &x, std::vector<unsigned> order) {
    assert(order.size() > 0);
    assert(x.size()%3 == 0);
    assert(order.size() <= x.size()/3);

    auto xnew = x;

    for (int i = 0; i < order.size()-1; ++i)
        xnew.segment(order[i]*3,3) = xnew.segment(order[i+1]*3,3);
    xnew.segment( order[order.size()-1]*3,3) = x.segment(order[0]*3,3);

    return xnew;
}


namespace PositionsVectorTransformer{

    struct AngleAxis{
    public:
        AngleAxis(double angle,Eigen::Vector3d axis)
                : angle_(angle), axis_(axis){};
        double angle_;
        Eigen::Vector3d axis_;
    };

    enum class RotationDirection {
        clockwise = 1,
        counterclockwise = -1
    };

    void translateCenterOfMassToOrigin(PositionsVector& positionsVector);

    void rotateAroundAxis(PositionsVector &positionsVector, double angle,
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
