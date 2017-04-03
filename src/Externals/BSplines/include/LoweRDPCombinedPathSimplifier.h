//
// Created by Michael Heuer on 29.08.16.
//

#include "LowePathSimplifier.h"
#include "RamerDouglasPeuckerPathSimplifier.h"

#ifndef RTQC_LOWERDPCOMBINEDSIMPLIFIER_H
#define RTQC_LOWERDPCOMBINEDSIMPLIFIER_H

/*! Convenience methods for Lowe and consecutive Ramer-Douglas-Peucker Simplification.
 * */
class LoweRDPCombinedPathSimplifier {
public:
  LoweRDPCombinedPathSimplifier();
  Eigen::MatrixXd simplifiedPath(const Eigen::MatrixXd &path,
                                 unsigned LoweMinSegLength,
                                 double RDPthreshold);

  Delib::MolecularTrajectory simplifiedPath(const Delib::MolecularTrajectory &path,
                                            unsigned LoweMinSegLength,
                                            double RDPThreshold);

  const Eigen::Matrix<unsigned,1,Eigen::Dynamic>& getSurvivingIndices(){ return survivingIndices_; };


private:
  LowePathSimplifier loweSimplifier_;
  RamerDouglasPeuckerPathSimplifier rdpSimplifier_;

  Eigen::Matrix<unsigned,1,Eigen::Dynamic> survivingIndices_;
  Eigen::MatrixXd simplifiedCoordMatrixAfterLoweAndRDP_;
};

#endif //RTQC_LOWERDPCOMBINEDSIMPLIFIER_H

