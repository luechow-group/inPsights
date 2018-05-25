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


    Eigen::MatrixXd correlationMatrix(const AllCentersSet& atomicA,
                                      const AllCentersSet& electronicA,
                                      const AllCentersSet& atomicB,
                                      const AllCentersSet& electronicB
    ){

        // TODO assert that A and B are subsets of particle pool

        auto M = ParticleKit::numberOfAtoms();
        auto N = ParticleKit::numberOfElectrons();

        auto MA = atomicA.size();
        auto MB = atomicB.size();

        auto NA = electronicA.size();
        auto NB = electronicB.size();

        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M+N,M+N);

        // Atoms A with Atoms B
        for (unsigned i = 0; i < MA; ++i)
            for (unsigned j = 0; j < MB; ++j)
                C(i,j) = LocalSimilarity::localSimilarity(atomicA[i], atomicB[j]);

        // Atoms A with Electrons B
        for (unsigned i = 0; i < MA; ++i)
            for (unsigned j = 0; j < NB; ++j)
                C(i,M+j) = LocalSimilarity::localSimilarity(atomicA[i], electronicB[j]);

        // Electrons A with Atoms B
        for (unsigned i = 0; i < NA; ++i)
            for (unsigned j = 0; j < MB; ++j)
                C(M+i,j) = LocalSimilarity::localSimilarity(electronicA[i],atomicB[j]);

        // Electrons A with Electrons B
        for (unsigned i = 0; i < NA; ++i)
            for (unsigned j = 0; j < NB; ++j)
                C(M+i,M+j) = LocalSimilarity::localSimilarity(electronicA[i],electronicB[j]);

        std::cout << std::setprecision(2) << C << std::endl<< std::endl;
        return C;
    }

    // compute only upper triangle and symmetrize
    static Eigen::MatrixXd correlationMatrix(const AllCentersSet& atomicA,
                                             const AllCentersSet& electronicA){

        auto M = ParticleKit::numberOfAtoms();
        auto N = ParticleKit::numberOfElectrons();

        auto MA = atomicA.size();
        auto NA = electronicA.size();


        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(M+N,M+N);

        // Atoms A with Atoms B
        for (unsigned i = 0; i < MA; ++i)
            for (unsigned j = i; j < MA; ++j)
                C(i,j) = LocalSimilarity::localSimilarity(atomicA[i], atomicA[j]);

        // Atoms A with Electrons A
        for (unsigned i = 0; i < MA; ++i)
            for (unsigned j = i; j < NA; ++j)
                C(i,M+j) = LocalSimilarity::localSimilarity(atomicA[i], electronicA[j]);

        // Electrons A with Atoms A
        for (unsigned i = 0; i < NA; ++i)
            for (unsigned j = i; j < MA; ++j)
                C(M+i,j) = LocalSimilarity::localSimilarity(electronicA[i],atomicA[j]);

        // Electrons A with Electrons A
        for (unsigned i = 0; i < NA; ++i)
            for (unsigned j = i; j < NA; ++j)
                C(M+i,M+j) = LocalSimilarity::localSimilarity(electronicA[i],electronicA[j]);

        std::cout << std::setprecision(2) << C << std::endl<< std::endl;
        // symmetrize the matrix
        for (unsigned i = 0; i < M+N; ++i)
            for (unsigned j = i+1; j < M+N; ++j)
                C(j,i) = C(i,j);

        return C;
    }


    double stucturalSimilarity(const MolecularGeometry& A,
                               const MolecularGeometry& B, double regularizationParameter) {

        NeighborhoodExpander expander;
        auto Aatomic = expander.computeExpansions(A,Type::Atomic);
        auto Aelectronic = expander.computeExpansions(A,Type::Electronic);

        auto Batomic = expander.computeExpansions(B,Type::Atomic);
        auto Belectronic = expander.computeExpansions(B,Type::Electronic);

        auto CAB = correlationMatrix(Aatomic,Aelectronic,Batomic,Belectronic);
        //auto bla = correlationMatrix(Aatomic,Aelectronic,Aatomic,Aelectronic);
        auto CAA = correlationMatrix(Aatomic,Aelectronic);
        auto CBB = correlationMatrix(Batomic,Belectronic);

        auto kAB = Sinkhorn::distance(CAB,regularizationParameter);
        auto kAA = Sinkhorn::distance(CAA,regularizationParameter);
        auto kBB = Sinkhorn::distance(CBB,regularizationParameter);

        return kAB/sqrt(kAA*kBB);
    }

};

#endif //AMOLQCPP_STRUCTURALSIMILARITY_H
