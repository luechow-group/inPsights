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

    Eigen::MatrixXd correlationMatrix(MolecularSpectrum& A, MolecularSpectrum& B) {
        auto N = ParticleKit::numberOfParticles();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);

        for (unsigned i = 0; i < N; ++i) {
            auto neA = ParticleKit::getNumberedTypeByIndex(i);
            if (!A.molecule_.findIndexByNumberedType(neA).first)
                continue;
            auto expA = A.molecularCenters_.find(neA)->second;

            for (unsigned j = 0; j < N; ++j) {
                auto neB = ParticleKit::getNumberedTypeByIndex(j);
                if (!B.molecule_.findIndexByNumberedType(neB).first)
                    continue;
                auto expB = B.molecularCenters_.find(neB)->second;

                C(i, j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }
        return C;
    }

    Eigen::MatrixXd correlationMatrixSame(MolecularSpectrum& A) {
        auto N = ParticleKit::numberOfParticles();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(N, N);

        for (unsigned i = 0; i < N; ++i) {
            auto neA = ParticleKit::getNumberedTypeByIndex(i);
            if (!A.molecule_.findIndexByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA)->second;

            for (unsigned j = i; j < N; ++j) {
                auto neB = ParticleKit::getNumberedTypeByIndex(j);
                if (!A.molecule_.findIndexByNumberedType(neB).first) continue;
                auto expB = A.molecularCenters_.find(neB)->second;

                auto val = LocalSimilarity::localSimilarity(expA, expB);
                C(i,j) = val;
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
        auto CAA = correlationMatrixSame(spectrumA);
        auto CBB = correlationMatrixSame(spectrumB);

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

};

#endif //AMOLQCPP_STRUCTURALSIMILARITY_H
