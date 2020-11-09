// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <omp.h>
#include "StructuralSimilarity.h"
#include "LocalSimilarity.h"
#include "Sinkhorn.h"

namespace SOAP {
    namespace StructuralSimilarity {

        Eigen::MatrixXd correlationMatrix(const MolecularSpectrum &A, const MolecularSpectrum &B) {
            assert(ParticleKit::isSubsetQ(A.molecule_)
                   && "The underlying molecule must be a subset of the selected particle kit.");
            assert(ParticleKit::isSubsetQ(B.molecule_)
                   && "The underlying molecule must be a subset of the selected particle kit.");

            auto N = ParticleKit::numberOfParticles();
            Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);

            EnumeratedType<int> enumType_i = {}, enumType_j = {};

#pragma omp parallel for default(none) shared(N, A, B, C, General::settings) private(enumType_i, enumType_j)
            for (unsigned i = 0; i < N; ++i) {
                //printf("Thread %d calculates correlation matrix elements\n", omp_get_thread_num());
                enumType_i = ParticleKit::getEnumeratedTypeByIndex(i);
                if (!A.molecule_.findIndexByEnumeratedType(enumType_i).first)
                    continue;
                const auto & expA = A.molecularCenters_.find(enumType_i)->second;

                for (unsigned j = 0; j < N; ++j) {
                    enumType_j = ParticleKit::getEnumeratedTypeByIndex(j);
                    if (!B.molecule_.findIndexByEnumeratedType(enumType_j).first)
                        continue;
                    const auto & expB = B.molecularCenters_.find(enumType_j)->second;
                    C(i, j) = LocalSimilarity::kernel(expA, expB, General::settings.zeta());
                }
            }
            return C;
        }

        Eigen::MatrixXd selfCorrelationMatrix(const MolecularSpectrum &A) {
            assert(ParticleKit::isSubsetQ(A.molecule_)
                   && "The underlying molecule must be a subset of the selected particle kit.");

            auto N = ParticleKit::numberOfParticles();
            Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);

            EnumeratedType<int> enumType_i = {}, enumType_j = {};

#pragma omp parallel for default(none) shared(N, A, C, General::settings) private(enumType_i, enumType_j)
            for (unsigned i = 0; i < N; ++i) {
                //printf("Thread %d calculates selfcorrelation matrix elements\n", omp_get_thread_num());
                enumType_i = ParticleKit::getEnumeratedTypeByIndex(i);
                if (!A.molecule_.findIndexByEnumeratedType(enumType_i).first)
                    continue;

                const auto& expA = A.molecularCenters_.find(enumType_i)->second;

                for (unsigned j = i; j < N; ++j) {
                    enumType_j = ParticleKit::getEnumeratedTypeByIndex(j);
                    if (!A.molecule_.findIndexByEnumeratedType(enumType_j).first)
                        continue;

                    const auto& expB = A.molecularCenters_.find(enumType_j)->second;

                    C(i, j) = LocalSimilarity::kernel(expA, expB, General::settings.zeta());
                }
            }
            // symmetrize the matrix
            C = C.selfadjointView<Eigen::Upper>();

            return C;
        }

        double kernel(const MolecularGeometry &A,
                      const MolecularGeometry &B, double gamma) {
            MolecularSpectrum spectrumA(A);
            MolecularSpectrum spectrumB(B);

            return kernel(spectrumA, spectrumB, gamma);
        }

        void spinFlippedSpectrum(const MolecularSpectrum &spectrumBspinflipped);

        double kernel(const MolecularSpectrum &spectrumA,
                      const MolecularSpectrum &spectrumB, double gamma) {

            auto CAB = correlationMatrix(spectrumA, spectrumB);
            auto CAA = selfCorrelationMatrix(spectrumA);
            auto CBB = selfCorrelationMatrix(spectrumB);

            auto kAB = Sinkhorn::distance(CAB, gamma);
            auto kAA = Sinkhorn::distance(CAA, gamma);
            auto kBB = Sinkhorn::distance(CBB, gamma);

            auto value = kAB / sqrt(kAA * kBB);

            if(General::settings.spinFlipCheck()) {
                auto spectrumBspinflipped = spectrumB;
                spectrumBspinflipped.flipSpins();

                // caluclate spin-flipped value
                auto CABsf = correlationMatrix(spectrumA, spectrumBspinflipped);
                auto kABsf  = Sinkhorn::distance(CABsf, gamma);
                auto valueSf = kABsf  / sqrt(kAA * kBB);
                if(value < valueSf)
                    value = valueSf;
            }

            return value;
        }



        double kernelDistance(const MolecularGeometry &A, const MolecularGeometry &B, double gamma) {
            return sqrt(2.0 - 2.0 * kernel(A, B, gamma));
        }

        double kernelDistance(const MolecularSpectrum &spectrumA, const MolecularSpectrum &spectrumB, double gamma) {
            return sqrt(2.0 - 2.0 * kernel(spectrumA, spectrumB, gamma));
        }
    }
}