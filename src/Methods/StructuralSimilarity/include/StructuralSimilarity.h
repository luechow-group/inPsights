//
// Created by Michael Heuer on 09.05.18.
//

#ifndef AMOLQCPP_STRUCTURALSIMILARITY_H
#define AMOLQCPP_STRUCTURALSIMILARITY_H

#include <Eigen/Core>
#include "ParticleKit.h"
#include "MolecularGeometry.h"
#include "LocalSimilarity.h"
#include "Sinkhorn.h"
#include "Environment.h"
#include <vector>
#include "NeighborhoodExpander.h"
#include "MolecularSpectrum.h"

namespace StructuralSimilarity{

    Eigen::MatrixXd correlationMatrix(MolecularSpectrum& A,
                                      MolecularSpectrum& B);

    Eigen::MatrixXd correlationMatrixSame(MolecularSpectrum& A);

    double stucturalSimilarity(const MolecularGeometry& A,
                               const MolecularGeometry& B, double regularizationParameter);
};

#endif //AMOLQCPP_STRUCTURALSIMILARITY_H
