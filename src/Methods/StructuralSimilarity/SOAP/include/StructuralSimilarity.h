// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_STRUCTURALSIMILARITY_H
#define INPSIGHTS_STRUCTURALSIMILARITY_H

#include <Eigen/Core>
#include "MolecularGeometry.h"
#include "MolecularSpectrum.h"
#include "SOAPSettings.h"

namespace SOAP {
    namespace StructuralSimilarity {

        Eigen::MatrixXd correlationMatrix(const MolecularSpectrum &A,
                                          const MolecularSpectrum &B);

        Eigen::MatrixXd selfCorrelationMatrix(const MolecularSpectrum &A);

        double kernel(const MolecularGeometry &A,
                      const MolecularGeometry &B, double gamma = General::settings.sinkhornGamma());

        double kernel(const MolecularSpectrum &spectrumA,
                      const MolecularSpectrum &spectrumB, double gamma = General::settings.sinkhornGamma());

        double kernelDistance(const MolecularGeometry &A,
                              const MolecularGeometry &B, double gamma = General::settings.sinkhornGamma());

        double kernelDistance(const MolecularSpectrum &spectrumA,
                              const MolecularSpectrum &spectrumB, double gamma = General::settings.sinkhornGamma());
    }
}

#endif //INPSIGHTS_STRUCTURALSIMILARITY_H
