//
// Created by Michael Heuer on 09.06.16.
//

#include "LowePathSimplifier.h"
#include "ContainerConverter.h"

LowePathSimplifier::LowePathSimplifier(){ }

std::vector<SegmentLowe> LowePathSimplifier::rateAndSplitSegment(SegmentLowe seg){

  unsigned idxOfPointWithMaximalSquaredDistanceToLine = 0;
  double squaredPointLineDistance=0;
  double maximalSquaredPointLineDistance=0;

  Eigen::VectorXd deltaStartEnd(path_.row(seg.endIdx_) - path_.row(seg.startIdx_));

  parametrizedLine_=Eigen::ParametrizedLine<double,Eigen::Dynamic>(path_.row(seg.startIdx_),(deltaStartEnd)
    .normalized());
  for (unsigned i = seg.startIdx_+1; i < seg.endIdx_-1; ++i) {
    squaredPointLineDistance = parametrizedLine_.squaredDistance(path_.row(i));
    if (squaredPointLineDistance >= maximalSquaredPointLineDistance){
      maximalSquaredPointLineDistance = squaredPointLineDistance;
      idxOfPointWithMaximalSquaredDistanceToLine = i;
    }
  }

  if (maximalSquaredPointLineDistance <= 1e-15) {
    if( deltaStartEnd.squaredNorm() <= 1e-15) seg.significance_=0;
    else seg.significance_ = std::numeric_limits<double>::infinity();

    return {seg};
  }
  else seg.significance_ = deltaStartEnd.squaredNorm()/maximalSquaredPointLineDistance;

  /* split this segment at the index with the maximal distance and deriveAndEvaluate */
  std::vector<SegmentLowe> subSegs,subSegs2;

  if ( (idxOfPointWithMaximalSquaredDistanceToLine-seg.startIdx_+1 >= minSegLength_) )
    subSegs = rateAndSplitSegment(SegmentLowe(seg.startIdx_, idxOfPointWithMaximalSquaredDistanceToLine));
  else
    subSegs = {SegmentLowe(seg.startIdx_,idxOfPointWithMaximalSquaredDistanceToLine,0)};

  if (seg.endIdx_-idxOfPointWithMaximalSquaredDistanceToLine+1 >= minSegLength_)
    subSegs2 = rateAndSplitSegment(SegmentLowe(idxOfPointWithMaximalSquaredDistanceToLine, seg.endIdx_));
  else
    subSegs2 ={SegmentLowe(idxOfPointWithMaximalSquaredDistanceToLine,seg.endIdx_,0)};


  /* concatenate vectors */
  subSegs.insert(subSegs.end(),
                std::make_move_iterator(subSegs2.begin()),
                std::make_move_iterator(subSegs2.end()));

  /* don't simplify when any sub-segment is more significant */
  for (auto &subSeg : subSegs){
    if (subSeg.significance_ >= seg.significance_) return subSegs;
  }
  return {seg};
}

Eigen::MatrixXd LowePathSimplifier::simplifiedPath(const Eigen::MatrixXd &path, unsigned minimumSegmentLength){
  minSegLength_ = minimumSegmentLength;

  // return the unchanged path if it is shorter than three segments
  if (path.rows() < 3) {
    survivingIndices_.resize(path.rows());
    for (unsigned i = 0; i < path.rows(); ++i) { survivingIndices_(i) = i; }
    return path;
  }

  assert(path.cols() >= 2 && "The path must be at least be two dimensional.");
  assert(minimumSegmentLength >= 2 && "The minimal segment length must be at least two.");

  path_ = path;


  SegmentLowe seg(0,unsigned(path_.rows())-1);
  std::vector<SegmentLowe> returnedSegments = rateAndSplitSegment(seg);

  /* construct simplified path from segment indices*/
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