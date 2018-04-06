//
// Created by Michael Heuer on 21.02.18.
//

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "ElectronicWaveFunctionProblem.h"
#include "CollectionParser.h"
#include "solver/newtonraphsonsolver.h"

int main(int argc, char *argv[]) {

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("BH3_Exp-em.wf");
    auto ac = electronicWaveFunctionProblem.getAtomsVector();
    std::cout << ac << std::endl;


    CollectionParser collectionParser;
    auto ec = collectionParser.electronsVectorFromJson("BH3_Max1.json");
    auto nsmooth = 2;
    auto x = ec.positionsVector().positionsAsEigenVector();

    std::cout << ec << std::endl;

    auto n = ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3;
    Eigen::VectorXd grad(n);
    electronicWaveFunctionProblem.putElectronsIntoNuclei(x,grad);

    std::cout << "gradient:" << std::endl;
    std::cout << ElectronsVector(grad,ec.spinTypesVector().spinTypesAsEigenVector()) << std::endl;

    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;
    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;


    Eigen::MatrixXd hess(n, n);
    electronicWaveFunctionProblem.hessian(x, hess);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(hess, true);
    std::cout << hess << std::endl;
    std::cout << eigenSolver.eigenvalues() << std::endl;
    std::cout << std::endl;
}
