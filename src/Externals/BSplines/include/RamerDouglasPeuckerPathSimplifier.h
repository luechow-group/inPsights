//
// Created by Michael Heuer on 01.06.16.
//

#ifndef RTQC_RAMERDOUGLASPEUCKERPATHSIMPLIFIER_H
#define RTQC_RAMERDOUGLASPEUCKERPATHSIMPLIFIER_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

namespace Delib{ class MolecularTrajectory; }

struct SegmentRDP {
  SegmentRDP(unsigned startIdx, unsigned endIdx)
    : startIdx_(startIdx),
      endIdx_(endIdx){};
  unsigned startIdx_,endIdx_;
};

/*! Implementation of the Ramer-Douglas-Peucker polyline simplification algorithm.
 * doi: 10.3138/FM57-6770-U75U-7727
 * */
class RamerDouglasPeuckerPathSimplifier {
public:

  RamerDouglasPeuckerPathSimplifier();

  Eigen::MatrixXd simplifiedPath(const Eigen::MatrixXd &path,const double threshold);
  Delib::MolecularTrajectory simplifiedPath(const Delib::MolecularTrajectory &path,const double threshold);

  const Eigen::Matrix<unsigned,1,Eigen::Dynamic>& getSurvivingIndices(){ return survivingIndices_; };

private:
  double threshold_;
  Eigen::ParametrizedLine<double,Eigen::Dynamic> parametrizedLine_;
  Eigen::MatrixXd path_,simplifiedPath_;
  Eigen::Matrix<unsigned,1,Eigen::Dynamic> survivingIndices_;
  std::vector<SegmentRDP> splitSegment(SegmentRDP seg);
};

#endif //RTQC_RAMERDOUGLASPEUCKERPATHSIMPLIFIER_H
