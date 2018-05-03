//
// Created by Michael Heuer on 02.05.18.
//

#ifndef AMOLQCPP_CUTOFF_H
#define AMOLQCPP_CUTOFF_H


#include <Eigen/Core>

class Cutoff{
public:
    explicit Cutoff(double cutoffRadius = 4., double cutoffWidth = 1., double centerWeight = 1.)
            : cutoffRadius_(cutoffRadius),
              cutoffWidth_(cutoffWidth), // GAP Paper: 1.0 Angstrom is regarded the length scale of atomic interactions
              centerWeight_(centerWeight),
              innerPlateauRadius_(cutoffRadius_ - cutoffWidth_)
    {}

    bool withinCutoffQ(const Eigen::Vector3d& position,
                       const Eigen::Vector3d& expansionCenter = Eigen::Vector3d::Zero()) const;

    
    double getWeight(double distanceFromExpansionCenter) const;

    double getWeight(const Eigen::Vector3d& position,
                     const Eigen::Vector3d& expansionCenter = Eigen::Vector3d::Zero()) const;

    Eigen::Vector3d getWeightGradient(const Eigen::Vector3d& position) const;


    double getCenterWeight(){ return centerWeight_; }

private:
    double distance(const Eigen::Vector3d &position,
                    const Eigen::Vector3d &expansionCenter = Eigen::Vector3d::Zero()) const;
    
    double cutoffRadius_, cutoffWidth_, centerWeight_, innerPlateauRadius_;
};

#endif //AMOLQCPP_CUTOFF_H
