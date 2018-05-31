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
        auto M = ParticleKit::numberOfAtoms();
        auto N = ParticleKit::numberOfElectrons();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M + N, M + N);

        // Atoms A with Atoms B
        for (unsigned i = 0; i < M; ++i) {

            auto neA = ParticleKit::getNumberedElementByIndex(i);
            if (!A.molecule_.atoms().findParticleByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA.toIntType())->second;

            for (unsigned j = 0; j < M; ++j) {

                auto neB = ParticleKit::getNumberedElementByIndex(j);
                if (!B.molecule_.atoms().findParticleByNumberedType(neB).first) continue;
                auto expB = B.molecularCenters_.find(neB.toIntType())->second;

                C(i, j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }

        // Atoms A with Electrons B
        for (unsigned i = 0; i < M; ++i) {

            auto neA = ParticleKit::getNumberedElementByIndex(i);
            if (!A.molecule_.atoms().findParticleByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA.toIntType())->second;

            for (unsigned j = 0; j < N; ++j) {
                auto neB = ParticleKit::getNumberedSpinByIndex(j);
                if (!B.molecule_.electrons().findParticleByNumberedType(neB).first) continue;
                auto expB = B.molecularCenters_.find(neB.toIntType())->second;

                C(i,M+j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }

        // Electrons A with Atoms B
        for (unsigned i = 0; i < N; ++i) {

            auto neA = ParticleKit::getNumberedSpinByIndex(i);
            if (!A.molecule_.electrons().findParticleByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA.toIntType())->second;

            for (unsigned j = 0; j < M; ++j) {
                auto neB = ParticleKit::getNumberedElementByIndex(j);
                if (!B.molecule_.atoms().findParticleByNumberedType(neB).first) continue;
                auto expB = B.molecularCenters_.find(neB.toIntType())->second;

                C(M + i, j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }

        // Electrons A with Electrons B
        for (unsigned i = 0; i < N; ++i) {

            auto neA = ParticleKit::getNumberedSpinByIndex(i);
            if (!A.molecule_.electrons().findParticleByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA.toIntType())->second;

            for (unsigned j = 0; j < N; ++j) {
                auto neB = ParticleKit::getNumberedSpinByIndex(j);
                if (!B.molecule_.electrons().findParticleByNumberedType(neB).first) continue;
                auto expB = B.molecularCenters_.find(neB.toIntType())->second;

                C(M+i,M+j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }
        return C;
    }

    Eigen::MatrixXd correlationMatrix(MolecularSpectrum& A) {
        auto M = ParticleKit::numberOfAtoms();
        auto N = ParticleKit::numberOfElectrons();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M + N, M + N);

        // Atoms A with Atoms A
        for (unsigned i = 0; i < M; ++i) {

            auto neA = ParticleKit::getNumberedElementByIndex(i);
            if (!A.molecule_.atoms().findParticleByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA.toIntType())->second;

            for (unsigned j = i; j < M; ++j) {

                auto neB = ParticleKit::getNumberedElementByIndex(j);
                if (!A.molecule_.atoms().findParticleByNumberedType(neB).first) continue;
                auto expB = A.molecularCenters_.find(neB.toIntType())->second;

                C(i, j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }

        // Atoms A with Electrons B
        for (unsigned i = 0; i < M; ++i) {

            auto neA = ParticleKit::getNumberedElementByIndex(i);
            if (!A.molecule_.atoms().findParticleByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA.toIntType())->second;

            for (unsigned j = 0; j < N; ++j) {
                auto neB = ParticleKit::getNumberedSpinByIndex(j);
                if (!A.molecule_.electrons().findParticleByNumberedType(neB).first) continue;
                auto expB = A.molecularCenters_.find(neB.toIntType())->second;

                C(i,M+j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }

        // Electrons A with Atoms B
        /*for (unsigned i = 0; i < N; ++i) {

            auto neA = ParticleKit::getNumberedSpinByIndex(i);
            if (!A.molecule_.electrons().findParticleByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA.toIntType())->second;

            for (unsigned j = 0; j < M; ++j) {
                auto neB = ParticleKit::getNumberedElementByIndex(j);
                if (!A.molecule_.atoms().findParticleByNumberedType(neB).first) continue;
                auto expB = A.molecularCenters_.find(neB.toIntType())->second;

                C(M+i,j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }*/

        // Electrons A with Electrons B
        for (unsigned i = 0; i < N; ++i) {

            auto neA = ParticleKit::getNumberedSpinByIndex(i);
            if (!A.molecule_.electrons().findParticleByNumberedType(neA).first) continue;
            auto expA = A.molecularCenters_.find(neA.toIntType())->second;

            for (unsigned j = i; j < N; ++j) {
                auto neB = ParticleKit::getNumberedSpinByIndex(j);
                if (!A.molecule_.electrons().findParticleByNumberedType(neB).first) continue;
                auto expB = A.molecularCenters_.find(neB.toIntType())->second;

                C(M+i,M+j) = LocalSimilarity::localSimilarity(expA, expB);
            }
        }

        // symmetrize the matrix
        for (unsigned i = 0; i < M+N; ++i)
            for (unsigned j = i+1; j < M+N; ++j)
                C(j,i) = C(i,j);

        return C;
    }

    double stucturalSimilarity(const MolecularGeometry& A,
                               const MolecularGeometry& B, double regularizationParameter) {

        MolecularSpectrum spectrumA(A);
        MolecularSpectrum spectrumB(B);

        auto CAB = correlationMatrix(spectrumA,spectrumB);
        //auto CAA = correlationMatrix(spectrumA,spectrumA);
        //auto CBB = correlationMatrix(spectrumB,spectrumB);
        auto CAA = correlationMatrix(spectrumA);
        auto CBB = correlationMatrix(spectrumB);

        //std::cout << CAB << std::endl << std::endl;
        //std::cout << CAA << std::endl << std::endl;
        //std::cout << CBB << std::endl << std::endl;

        auto kAB = Sinkhorn::distance(CAB,regularizationParameter);
        auto kAA = Sinkhorn::distance(CAA,regularizationParameter);
        auto kBB = Sinkhorn::distance(CBB,regularizationParameter);

        //std::cout << kAB << std::endl;
        //std::cout << kAA << std::endl;
        //std::cout << kBB << std::endl;

        return kAB/sqrt(kAA*kBB);
    }

};

#endif //AMOLQCPP_STRUCTURALSIMILARITY_H
