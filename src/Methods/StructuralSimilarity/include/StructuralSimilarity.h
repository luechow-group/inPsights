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

namespace StructuralSimilarity{


    static Eigen::MatrixXd correlationMatrix(const MolecularGeometry& A, const MolecularGeometry& B){

        // TODO assert that A and B are subsets of particle pool

        auto M = ParticleKit::numberOfAtoms();
        auto N = ParticleKit::numberOfElectrons();

        auto MA = A.atoms().numberOfEntities();
        auto MB = B.atoms().numberOfEntities();

        auto NA = A.electrons().numberOfEntities();
        auto NB = B.electrons().numberOfEntities();

        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M+N,M+N);

        // TODO Do LocalSimilarity::computeExpansion here once
        //NeighborhoodExpander neighborhoodExpander;
        //std::vector<std::map<int, NeighborhoodExpansion>> expansionsA;
        //for (int k = 0; k < NA; ++k) {
        //    expansionsA.emplace(N); // Exp for each particle, - what for doubled particles?
        //}

        // Environments
        // Atoms A with Atoms B
        for (unsigned i = 0; i < MA; ++i) {
            for (unsigned j = 0; j < MB; ++j) {
                Environment eAi(A,A.atoms()[i].position());
                Environment eBj(B,A.atoms()[j].position());
                C(i,j) = LocalSimilarity::localSimilarity(eAi,eBj);
            }
        }

        // Atoms A with Electrons B
        for (unsigned i = 0; i < MA; ++i) {
            for (unsigned j = 0; j < NB; ++j) {
                Environment eAi(A,A.atoms()[i].position());
                Environment eBj(B,A.electrons()[j].position());
                C(i,M+j) = LocalSimilarity::localSimilarity(eAi,eBj);
            }
        }

        // Electrons A with Atoms B
        for (unsigned i = 0; i < NA; ++i) {
            for (unsigned j = 0; j < MB; ++j) {
                Environment eAi(A,A.electrons()[i].position());
                Environment eBj(B,A.atoms()[j].position());
                C(M+i,j) = LocalSimilarity::localSimilarity(eAi,eBj);
            }
        }

        // Electrons A with Electrons B
        for (unsigned i = 0; i < NA; ++i) {
            for (unsigned j = 0; j < NB; ++j) {
                Environment eAi(A,A.electrons()[i].position());
                Environment eBj(B,A.electrons()[j].position());
                C(M+i,M+j) = LocalSimilarity::localSimilarity(eAi,eBj);
            }
        }
        return C;
    }

    // compute only upper triangle and symmetrize
    static Eigen::MatrixXd correlationMatrix(const MolecularGeometry& A){


        auto M = ParticleKit::numberOfAtoms();
        auto N = ParticleKit::numberOfElectrons();

        auto MA = A.atoms().numberOfEntities();

        auto NA = A.electrons().numberOfEntities();

        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M+N,M+N);


        // Environments
        // Atoms A with Atoms B
        for (unsigned i = 0; i < MA; ++i) {
            for (unsigned j = i; j < MA; ++j) {
                Environment eAi(A,A.atoms()[i].position());
                Environment eAj(A,A.atoms()[j].position());
                C(i,j) = LocalSimilarity::localSimilarity(eAi,eAj);
            }
        }

        // Atoms A with Electrons B
        for (unsigned i = 0; i < MA; ++i) {
            for (unsigned j = i; j < NA; ++j) {
                Environment eAi(A,A.atoms()[i].position());
                Environment eAj(A,A.electrons()[j].position());
                C(i,M+j) = LocalSimilarity::localSimilarity(eAi,eAj);
            }
        }

        // Electrons A with Atoms B
        for (unsigned i = 0; i < NA; ++i) {
            for (unsigned j = i; j < MA; ++j) {
                Environment eAi(A,A.electrons()[i].position());
                Environment eAj(A,A.atoms()[j].position());
                C(M+i,j) = LocalSimilarity::localSimilarity(eAi,eAj);
            }
        }

        // Electrons A with Electrons B
        for (unsigned i = 0; i < NA; ++i) {
            for (unsigned j = i; j < NA; ++j) {
                Environment eAi(A,A.electrons()[i].position());
                Environment eAj(A,A.electrons()[j].position());
                C(M+i,M+j) = LocalSimilarity::localSimilarity(eAi,eAj);
            }
        }

        // symmetrize the matrix
        for (unsigned i = 0; i < M; ++i) {
            for (unsigned j = i+1; j < M; ++j) {
                C(j,i) = C(i,j);
            }
        }

        return C;
    }


    static double stucturalSimilarity(const MolecularGeometry& A,
                                      const MolecularGeometry& B, double regularizationParameter) {

        auto kAB = Sinkhorn::distance(correlationMatrix(A,B),regularizationParameter);
        auto kAA = Sinkhorn::distance(correlationMatrix(A),regularizationParameter);
        auto kBB = Sinkhorn::distance(correlationMatrix(B),regularizationParameter);

        return kAB/sqrt(kAA*kBB);
    }

};

#endif //AMOLQCPP_STRUCTURALSIMILARITY_H
