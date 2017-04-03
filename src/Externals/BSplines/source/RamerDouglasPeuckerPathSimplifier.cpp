//
// Created by Michael Heuer on 01.06.16.
//

#include "RamerDouglasPeuckerPathSimplifier.h"
#include "ContainerConverter.h"

RamerDouglasPeuckerPathSimplifier::RamerDouglasPeuckerPathSimplifier(){ }

std::vector<SegmentRDP> RamerDouglasPeuckerPathSimplifier::splitSegment(SegmentRDP seg) {
  unsigned dMaxIdx = 0;
  double d=0;
  double dMax = 0;

  Eigen::VectorXd deltaStartEnd(path_.row(seg.endIdx_) -path_.row(seg.startIdx_));

  parametrizedLine_ = Eigen::ParametrizedLine<double,Eigen::Dynamic>(path_.row(seg.startIdx_),(deltaStartEnd)
    .normalized());
  for (unsigned i = seg.startIdx_+1; i < seg.endIdx_; ++i) {
    d = parametrizedLine_.distance(path_.row(i));
    if (d >= dMax){
      dMax = d;
      dMaxIdx = i;
    }
  }

  if (dMax > threshold_) {
    std::vector<SegmentRDP> result1,result2;
    result1  = splitSegment(SegmentRDP(seg.startIdx_,dMaxIdx));
    result2 = splitSegment(SegmentRDP(dMaxIdx,seg.endIdx_));
    result1.insert(result1.end(),
                  std::make_move_iterator(result2.begin()),
                  std::make_move_iterator(result2.end()));
    return result1;
  }
  else return {SegmentRDP(seg.startIdx_,seg.endIdx_)};
}


Eigen::MatrixXd
RamerDouglasPeuckerPathSimplifier::simplifiedPath(const Eigen::MatrixXd &path, const double threshold) {
  threshold_ = std::abs(threshold);

  // return the unchanged path if it is shorter than three segments
  if (path.rows() < 3) {
    survivingIndices_.resize(path.rows());
    for (unsigned i = 0; i < path.rows(); ++i) { survivingIndices_(i) = i; }
    return path;
  }

  assert(path.cols() >= 2 && "The path must be at least be two dimensional.");

  path_ = path;

  SegmentRDP seg(0,unsigned(path_.rows())-1);
  std::vector<SegmentRDP>returnedSegments = splitSegment(seg);

  /* construct simplified path from egment indices*/
  survivingIndices_.resize(returnedSegments.size()+1);
  simplifiedPath_.resize(returnedSegments.size()+1,path_.cols());

  for (unsigned i = 0; i < returnedSegments.size(); ++i) {
    survivingIndices_(i) = returnedSegments[i].startIdx_;
    simplifiedPath_.row(i) = path_.row(returnedSegments[i].startIdx_);
  }

  survivingIndices_.tail(1)(0) = returnedSegments.back().endIdx_;
  simplifiedPath_.bottomRows(1).row(0) = path_.row(returnedSegments.back().endIdx_);

  return simplifiedPath_;
}
