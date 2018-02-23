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
    auto ac = electronicWaveFunctionProblem.getAtomCollection();
    std::cout << ac << std::endl;


    CollectionParser collectionParser;
    auto ec = collectionParser.electronCollectionFromJson("BH3_Max1.json");
    auto nsmooth = 2;
    auto x = ec.positionsAsEigenVector();


    auto n = ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3;
    Eigen::VectorXd grad(n);
    electronicWaveFunctionProblem.putElectronsIntoNuclei(x,grad);

    std::cout << ElectronCollection(grad,ec.spinTypesAsEigenVector()) << std::endl;

    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;
    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;

    std::vector<double >epsilons({1E-2,1E-3,1E-4,1E-5,1E-6,1E-7,1E-8,1E-9,1E-10,1E-11});

    auto eps= std::numeric_limits<double >::epsilon();

    std::cout << eps << std::endl;
    for (int j = 0; j < x.size(); ++j) {
        std::cout << "x_i = " << x(j) << ", h_i = " << std::sqrt(eps) * (std::abs(x(j)) + std::sqrt(eps)) << std::endl;
    }

    std::cout << "default eps = " << electronicWaveFunctionProblem.gethDelta() << " " << x(1) << std::endl;

    for (auto & i : epsilons) {
        electronicWaveFunctionProblem.sethDelta(i);
        std::cout << "eps = " << electronicWaveFunctionProblem.gethDelta() << std::endl;
        Eigen::MatrixXd hess(n, n);
        electronicWaveFunctionProblem.hessian(x, hess);
        //std::cout << hess << std::endl;

        auto relevantBlock = hess.block((8 - nsmooth) * 3, (8 - nsmooth) * 3, 3 * nsmooth, 3 * nsmooth);
        Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(relevantBlock, true);
        std::cout << relevantBlock << std::endl;
        //Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(hess,true);
        auto eigenvalues = eigenSolver.eigenvalues();
        std::cout << eigenvalues << std::endl;
        //auto eigenvectors = eigenSolver.eigenvectors();
        //std::cout << eigenvectors << std::endl;
        std::cout << std::endl;
    }
    /*
    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
    crit.iterations = 1000;
    crit.gradNorm = 1e-8;
    cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);
    Eigen::VectorXd x = x0;
    solver.minimize(electronicWaveFunctionProblem, x);
    std::cout << "max: " << x << std::endl;
     */

    //cppoptlib::NewtonRaphsonSolver<ElectronicWaveFunctionProblem> solver;
    //solver.setDebug(cppoptlib::DebugLevel::High);
    //cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
    //crit.iterations = 20;
    //crit.gradNorm = 1e-6;
    //solver.setStopCriteria(crit);
    //solver.minimize(f,x);


    //CollectionParser collectionParser;
    ////auto ecA = collectionParser.electronCollectionFromJson("Ethane-glob-max.json");
    //auto ecA = ElectronCollection(x0,Eigen::Vector2i(1,-1));
    //auto ecB = ecA;
    //std::cout << ecA << std::endl;



}
