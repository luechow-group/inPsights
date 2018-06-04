//
// Created by Michael Heuer on 02.06.18.
//

#include "StructuralSimilarity.h"

namespace StructuralSimilarity{

    Eigen::MatrixXd correlationMatrix(MolecularSpectrum& A, MolecularSpectrum& B) {
        auto N = ParticleKit::numberOfParticles();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);

        for (unsigned i = 0; i < N; ++i) {
            auto numberedType_i = ParticleKit::getNumberedTypeByIndex(i);
            if (!A.molecule_.findIndexByNumberedType(numberedType_i).first)
                continue;
            auto expA = A.molecularCenters_.find(numberedType_i)->second;

            for (unsigned j = 0; j < N; ++j) {
                auto numberedType_j = ParticleKit::getNumberedTypeByIndex(j);
                if (!B.molecule_.findIndexByNumberedType(numberedType_j).first)
                    continue;
                auto expB = B.molecularCenters_.find(numberedType_j)->second;

                C(i, j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }
        return C;
    }

    Eigen::MatrixXd selfCorrelationMatrix(MolecularSpectrum &A) {
        auto N = ParticleKit::numberOfParticles();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);

        for (unsigned i = 0; i < N; ++i) {
            auto numberedType_i = ParticleKit::getNumberedTypeByIndex(i);
            if (!A.molecule_.findIndexByNumberedType(numberedType_i).first) continue;
            auto expA = A.molecularCenters_.find(numberedType_i)->second;

            for (unsigned j = i; j < N; ++j) {
                auto numberedType_j = ParticleKit::getNumberedTypeByIndex(j);
                if (!A.molecule_.findIndexByNumberedType(numberedType_j).first) continue;
                auto expB = A.molecularCenters_.find(numberedType_j)->second;

                C(i,j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }
        // symmetrize the matrix
        for (unsigned i = 0; i < N; ++i)
            for (unsigned j = i+1; j < N; ++j)
                C(j,i) = C(i,j);

        return C;
    }

    double stucturalSimilarity(const MolecularGeometry& A,
                               const MolecularGeometry& B, double regularizationParameter) {

        MolecularSpectrum spectrumA(A);
        MolecularSpectrum spectrumB(B);

        auto CAB = correlationMatrix(spectrumA,spectrumB);
        auto CAA = selfCorrelationMatrix(spectrumA);
        auto CBB = selfCorrelationMatrix(spectrumB);

        std::cout << CAB << std::endl << std::endl;
        std::cout << CAA << std::endl << std::endl;
        std::cout << CBB << std::endl << std::endl;

        auto kAB = Sinkhorn::distance(CAB,regularizationParameter);
        auto kAA = Sinkhorn::distance(CAA,regularizationParameter);
        auto kBB = Sinkhorn::distance(CBB,regularizationParameter);

        std::cout << kAB << std::endl;
        std::cout << kAA << std::endl;
        std::cout << kBB << std::endl;

        return kAB/sqrt(kAA*kBB);
    }

}