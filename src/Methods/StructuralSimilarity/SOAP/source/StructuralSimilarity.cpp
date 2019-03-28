//
// Created by Michael Heuer on 02.06.18.
//
#include <omp.h>
#include "StructuralSimilarity.h"
#include "LocalSimilarity.h"
#include "Sinkhorn.h"

namespace StructuralSimilarity{

    Eigen::MatrixXd correlationMatrix(const MolecularSpectrum& A, const MolecularSpectrum& B) {
        assert(ParticleKit::isSubsetQ(A.molecule_)
               && "The underlying molecule must be a subset of the selected particle kit.");
        assert(ParticleKit::isSubsetQ(B.molecule_)
               && "The underlying molecule must be a subset of the selected particle kit.");

        auto N = ParticleKit::numberOfParticles();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);
        NumberedType<int> numberedType_i, numberedType_j;
        TypeSpecificNeighborhoodsAtOneCenter expA,expB;

        #pragma omp parallel for default(none) shared(N,A,B,C,SOAPExpansion::settings) private(numberedType_i,numberedType_j,expA,expB)
        for (unsigned i = 0; i < N; ++i) {
            //printf("Thread %d calculates correlation matrix elements\n", omp_get_thread_num());
            numberedType_i = ParticleKit::getNumberedTypeByIndex(i);
            if (!A.molecule_.findIndexByNumberedType(numberedType_i).first)
                continue;
            expA = A.molecularCenters_.find(numberedType_i)->second;

            for (unsigned j = 0; j < N; ++j) {
                numberedType_j = ParticleKit::getNumberedTypeByIndex(j);
                if (!B.molecule_.findIndexByNumberedType(numberedType_j).first)
                    continue;
                expB = B.molecularCenters_.find(numberedType_j)->second;
                C(i, j) = LocalSimilarity::kernel(expA, expB, SOAPExpansion::settings.zeta());
            }
        }
        return C;
    }

    Eigen::MatrixXd selfCorrelationMatrix(const MolecularSpectrum &A) {
        assert(ParticleKit::isSubsetQ(A.molecule_)
               && "The underlying molecule must be a subset of the selected particle kit.");

        auto N = ParticleKit::numberOfParticles();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);

        NumberedType<int> numberedType_i, numberedType_j;
        TypeSpecificNeighborhoodsAtOneCenter expA,expB;
        #pragma omp parallel for default(none) shared(N,A,C, SOAPExpansion::settings) private(numberedType_i,numberedType_j,expA,expB)
        for (unsigned i = 0; i < N; ++i) {
            //printf("Thread %d calculates selfcorrelation matrix elements\n", omp_get_thread_num());
            numberedType_i = ParticleKit::getNumberedTypeByIndex(i);
            if (!A.molecule_.findIndexByNumberedType(numberedType_i).first) continue;
            expA = A.molecularCenters_.find(numberedType_i)->second;

            for (unsigned j = i; j < N; ++j) {
                numberedType_j = ParticleKit::getNumberedTypeByIndex(j);
                if (!A.molecule_.findIndexByNumberedType(numberedType_j).first) continue;
                expB = A.molecularCenters_.find(numberedType_j)->second;

                C(i,j) = LocalSimilarity::kernel(expA, expB, SOAPExpansion::settings.zeta());
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