/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
