//
// Created by Michael Heuer on 29.08.16.
//

#include "LoweRDPCombinedPathSimplifier.h"
#include "ContainerConverter.h"

LoweRDPCombinedPathSimplifier::LoweRDPCombinedPathSimplifier()
  :loweSimplifier_(), rdpSimplifier_() {}

Eigen::MatrixXd
LoweRDPCombinedPathSimplifier::simplifiedPath(const Eigen::MatrixXd &path, unsigned int LoweMinSegLength,
                                              double RDPThreshold) {

  auto simplifiedCoordMatrixAfterLowe = loweSimplifier_.simplifiedPath(path,LoweMinSegLength);
  auto remainingIndicesFromOriginalPathAfterLowe = loweSimplifier_.getSurvivingIndices();

  simplifiedCoordMatrixAfterLoweAndRDP_ = rdpSimplifier_.simplifiedPath(simplifiedCoordMatrixAfterLowe,RDPThreshold);
  auto remainingIndicesFromLoweSimplifiedPathAfterRDP = rdpSimplifier_.getSurvivingIndices();

  Eigen::Matrix<unsigned,1,Eigen::Dynamic>
    remainingIndicesFromOriginalPathAfterLoweAndRDP(simplifiedCoordMatrixAfterLoweAndRDP_.rows());

  // pick indices in original path
  assert(remainingIndicesFromLoweSimplifiedPathAfterRDP.cols() == simplifiedCoordMatrixAfterLoweAndRDP_.rows());
  for (int i = 0; i <remainingIndicesFromLoweSimplifiedPathAfterRDP.cols(); ++i) {
    auto indexInLoweSimplifedPath = remainingIndicesFromLoweSimplifiedPathAfterRDP[i];
    remainingIndicesFromOriginalPathAfterLoweAndRDP[i] = remainingIndicesFromOriginalPathAfterLowe[indexInLoweSimplifedPath];
  }

  survivingIndices_ = remainingIndicesFromLoweSimplifiedPathAfterRDP;
  return simplifiedCoordMatrixAfterLoweAndRDP_;
}
