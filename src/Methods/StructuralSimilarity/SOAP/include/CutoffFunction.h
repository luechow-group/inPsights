//
// Created by Michael Heuer on 02.05.18.
//

#ifndef INPSIGHTS_CUTOFF_H
#define INPSIGHTS_CUTOFF_H


#include <Eigen/Core>

//TODO make Cutoff a namespace and put cutoff settings into Expansion settings class
namespace CutoffFunction{
    bool withinCutoffRadiusQ(double distance);

    double getWeight(double distanceFromExpansionCenter);

    double getWeight(const Eigen::Vector3d& position,
                     const Eigen::Vector3d& expansionCenter = Eigen::Vector3d::Zero());

    Eigen::Vector3d getWeightGradient(const Eigen::Vector3d& position);

    double getCenterWeight();

    double distance(const Eigen::Vector3d &position,
                    const Eigen::Vector3d &expansionCenter = Eigen::Vector3d::Zero());
};

#endif //INPSIGHTS_CUTOFF_H
