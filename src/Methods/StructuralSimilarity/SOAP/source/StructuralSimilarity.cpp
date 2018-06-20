//
// Created by Michael Heuer on 02.06.18.
//

#include "StructuralSimilarity.h"

namespace StructuralSimilarity{

    Eigen::MatrixXd correlationMatrix(const MolecularSpectrum& A, const MolecularSpectrum& B) {
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

                C(i, j) = LocalSimilarity::kernel(expA, expB);
            }
        }
        return C;
    }

    Eigen::MatrixXd selfCorrelationMatrix(const MolecularSpectrum &A) {
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

                C(i,j) = LocalSimilarity::kernel(expA, expB);
            }
        }
        // symmetrize the matrix
        for (unsigned i = 0; i < N; ++i)
            for (unsigned j = i+1; j < N; ++j)
                C(j,i) = C(i,j);

        return C;
    }

    double kernel(const MolecularGeometry &A,
                  const MolecularGeometry &B, double gamma) {
        MolecularSpectrum spectrumA(A);
        MolecularSpectrum spectrumB(B);

        return kernel(spectrumA, spectrumB, gamma);
    }

    double kernel(const MolecularSpectrum &spectrumA,
                  const MolecularSpectrum &spectrumB, double gamma) {

        auto CAB = correlationMatrix(spectrumA,spectrumB);
        auto CAA = selfCorrelationMatrix(spectrumA);
        auto CBB = selfCorrelationMatrix(spectrumB);

        //std::cout << CAB << std::endl << std::endl;
        //std::cout << CAA << std::endl << std::endl;
        //std::cout << CBB << std::endl << std::endl;


        double eps = 1e-8; //TODO move in settings?
        auto kAB = Sinkhorn::distance(CAB,gamma,eps);
        auto kAA = Sinkhorn::distance(CAA,gamma,eps);
        auto kBB = Sinkhorn::distance(CBB,gamma,eps);

        //std::cout << kAB << std::endl;
        //std::cout << kAA << std::endl;
        //std::cout << kBB << std::endl;

        return kAB/sqrt(kAA*kBB);
    }

    double kernelDistance(const MolecularGeometry &A, const MolecularGeometry &B, double gamma) {
        return sqrt(2.0-2.0*kernel(A, B, gamma));
    }

    double kernelDistance(const MolecularSpectrum& spectrumA, const MolecularSpectrum& spectrumB, double gamma) {
        return sqrt(2.0-2.0*kernel(spectrumA, spectrumB, gamma));
    }
}