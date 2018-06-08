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

    Eigen::MatrixXd correlationMatrix(const MolecularSpectrum& A,
                                      const MolecularSpectrum& B);

    Eigen::MatrixXd selfCorrelationMatrix(const MolecularSpectrum &A);

    double kernel(const MolecularGeometry &A,
                  const MolecularGeometry &B, double gamma = ExpansionSettings::gamma);

    double kernel(const MolecularSpectrum &spectrumA,
                  const MolecularSpectrum &spectrumB, double gamma = ExpansionSettings::gamma);

    double kernelDistance(const MolecularGeometry &A,
                          const MolecularGeometry &B, double gamma = ExpansionSettings::gamma);
};

#endif //AMOLQCPP_STRUCTURALSIMILARITY_H
