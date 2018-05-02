//
// Created by Michael Heuer on 02.05.18.
//

#ifndef AMOLQCPP_CUTOFF_H
#define AMOLQCPP_CUTOFF_H

#include <cmath>
#include <Eigen/Core>

class Cutoff{
public:
    explicit Cutoff(double rCut = 4., double rCutWidth = 1./2)
            : rCut_(rCut),
              rCutWidth_(rCutWidth)
    {}

    double getWeight(double r){
        double innerPlateauRadius = (rCut_-rCutWidth_);

        if (0 < r && r <= innerPlateauRadius)
            return 1.;
        else if (innerPlateauRadius < r && r <= rCut_)
            return 0.5*( 1 + cos( M_PI*(r-innerPlateauRadius)/rCutWidth_) );
        else
            return 0.;
    };

    double getWeight(const Eigen::Vector3d& rVec, const Eigen::Vector3d& rCenter = Eigen::Vector3d::Zero()){
        return getWeight((rVec-rCenter).eval().norm());
    }

    Eigen::Vector3d getWeightGradient(const Eigen::Vector3d& rVec){

        double r = rVec.norm();
        Eigen::Vector3d rDir = rVec.normalized();

        double innerPlateauRadius = (rCut_-rCutWidth_);

        if (r <= innerPlateauRadius || r > rCut_)
            return Eigen::Vector3d::Zero();
        else if (innerPlateauRadius < r && r <= rCut_)
            return 0.5*( 1 + sin( M_PI*(r-innerPlateauRadius)/rCutWidth_)*M_PI/rCutWidth_ )*rDir;
    };


private:
    double rCut_, rCutWidth_;
};

#endif //AMOLQCPP_CUTOFF_H
