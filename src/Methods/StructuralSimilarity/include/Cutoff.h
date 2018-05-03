//
// Created by Michael Heuer on 02.05.18.
//

#ifndef AMOLQCPP_CUTOFF_H
#define AMOLQCPP_CUTOFF_H


#include <Eigen/Core>

class Cutoff{
public:
    explicit Cutoff(double cutoffRadius = 4., double cutoffWidth = 1., double centerWeight = 1.);

    bool withinCutoffRadiusQ(double distance) const;

    
    double getWeight(double distanceFromExpansionCenter) const;

    double getWeight(const Eigen::Vector3d& position,
                     const Eigen::Vector3d& expansionCenter = Eigen::Vector3d::Zero()) const;

    Eigen::Vector3d getWeightGradient(const Eigen::Vector3d& position) const;

    double getCenterWeight();

    static double distance(const Eigen::Vector3d &position,
                           const Eigen::Vector3d &expansionCenter = Eigen::Vector3d::Zero());
private:
    double cutoffRadius_, cutoffWidth_, centerWeight_, innerPlateauRadius_;
};

#endif //AMOLQCPP_CUTOFF_H
