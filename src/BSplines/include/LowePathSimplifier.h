//
// Created by Michael Heuer on 09.06.16.
//

#ifndef RTQC_LOWEPATHSIMPLIFIER_H
#define RTQC_LOWEPATHSIMPLIFIER_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <limits>
#include <vector>


struct SegmentLowe {
  SegmentLowe(unsigned first, unsigned last, double significance = 0)
    : startIdx_(first),
      endIdx_(last),
      significance_(significance) {};
  unsigned startIdx_,endIdx_;
  double significance_;
};

/*! Implementation of the Lowe algorithm. The algorithm simplifies a path of any dimensionality by a significance criterion
 * employing the distance between the segment ends and the maximum deviation of the line connecting them.
 * doi: 10.1016/0004-3702(87)90070-1
 * Deviating from the original algorithm we do not employ the last step where all segments with significances <= 4 are
 * deleted because we do not find it necessary for our application and it would involve more computational effort1k.
 * */
class LowePathSimplifier {
public:
  LowePathSimplifier();

  Eigen::MatrixXd simplifiedPath(const Eigen::MatrixXd &path, unsigned minimumSegmentLength = 3);

  const Eigen::Matrix<unsigned,1,Eigen::Dynamic>& getSurvivingIndices(){ return survivingIndices_; };

private:
  unsigned minSegLength_;
  Eigen::ParametrizedLine<double,Eigen::Dynamic> parametrizedLine_;
  Eigen::MatrixXd path_,simplifiedPath_;
  Eigen::Matrix<unsigned,1,Eigen::Dynamic> survivingIndices_;

  std::vector<SegmentLowe> rateAndSplitSegment(SegmentLowe seg);
};

#endif //RTQC_LOWEPATHSIMPLIFIER_H
