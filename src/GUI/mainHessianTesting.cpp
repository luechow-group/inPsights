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
    auto ec = collectionParser.electronCollectionFromJson("BH3_Max2.json");
    auto x = ec.positionsAsEigenVector();

    //ElectronCollection ec;
    //ec.append(Electron(ac[0],Spin::SpinType::alpha));
    //ec.append(Electron(ac[1],Spin::SpinType::alpha));
    //ec.append(Electron(ac[2],Spin::SpinType::alpha));
    //ec.append(Electron(ac[3],Spin::SpinType::alpha));
    //ec.append(Electron(ac[0],Spin::SpinType::beta));
    //ec.append(Electron(ac[1].position().array()/2.0, Spin::SpinType::beta));
    //ec.append(Electron(ac[2].position().array()/2.0, Spin::SpinType::beta));
    //ec.append(Electron(ac[3].position().array()/2.0, Spin::SpinType::beta));
    //auto x = ec.positionsAsEigenVector();
    //std::cout << ec << std::endl;

    auto n = ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3;
    Eigen::VectorXd grad(n);
    electronicWaveFunctionProblem.putElectronsIntoNuclei(x,grad);
    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;
    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;


    Eigen::MatrixXd hess(n,n);
    electronicWaveFunctionProblem.hessian(x,hess);
    std::cout << hess << std::endl;

    Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(hess, false);
    //Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(hess,true);
    auto eigenvalues = eigenSolver.eigenvalues();
    std::cout << eigenvalues << std::endl;
    //auto eigenvectors = eigenSolver.eigenvectors();
    //std::cout << eigenvectors << std::endl;

    //Eigen::VectorXd x(2*3);
    //x <<
    //   0,0,-0.700144,\
    //   0,0,+0.700144;

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
