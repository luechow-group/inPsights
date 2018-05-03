//
// Created by Michael Heuer on 03.05.18.
//

#include "Cutoff.h"
#include <cmath>

Cutoff::Cutoff(double cutoffRadius, double cutoffWidth, double centerWeight)
        : cutoffRadius_(cutoffRadius),
          cutoffWidth_(cutoffWidth),
          centerWeight_(centerWeight),
          innerPlateauRadius_(cutoffRadius_ - cutoffWidth_)
{}

bool Cutoff::withinCutoffRadiusQ(double distance) const {
    return distance < cutoffRadius_;
}

double Cutoff::getWeight(double distanceFromExpansionCenter) const{
    //TODO delete centerWeight and use: 'if (0 <= distanceFromExpansionCenter...' instead?
    if (0 < distanceFromExpansionCenter && distanceFromExpansionCenter <= innerPlateauRadius_)
        return 1.;
    else if (innerPlateauRadius_ < distanceFromExpansionCenter && distanceFromExpansionCenter <= cutoffRadius_)
        return 0.5*( 1 + cos( M_PI*(distanceFromExpansionCenter-innerPlateauRadius_)/cutoffWidth_) );
    else
        return 0.;
};


double Cutoff::getWeight(const Eigen::Vector3d& position,
                         const Eigen::Vector3d& expansionCenter) const {
    return getWeight(distance(position, expansionCenter));
}

Eigen::Vector3d Cutoff::getWeightGradient(const Eigen::Vector3d&position ) const {

    double distanceFromExpansionCenter =position .norm();
    Eigen::Vector3d direction =position .normalized();

    if (distanceFromExpansionCenter <= innerPlateauRadius_ || distanceFromExpansionCenter > cutoffRadius_)
        return Eigen::Vector3d::Zero();
    else if (innerPlateauRadius_ < distanceFromExpansionCenter && distanceFromExpansionCenter <= cutoffRadius_)
        return 0.5*( 1 + sin( M_PI*(distanceFromExpansionCenter-innerPlateauRadius_)/cutoffWidth_)*M_PI/cutoffWidth_ )*direction;
};


double Cutoff::getCenterWeight(){ return centerWeight_; }


double Cutoff::distance(const Eigen::Vector3d &position,
                        const Eigen::Vector3d &expansionCenter){
    return (position-expansionCenter).eval().norm();
}
