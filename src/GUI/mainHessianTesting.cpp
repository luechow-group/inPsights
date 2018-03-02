//
// Created by Michael Heuer on 21.02.18.
//

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <QApplication>

#include "MoleculeWidget.h"
#include "ElectronCollection3D.h"
#include "AtomCollection3D.h"
#include "Polyline.h"
#include "ElectronicWaveFunctionProblem.h"
#include "CollectionParser.h"
#include "solver/newtonraphsonsolver.h"


int main(int argc, char *argv[]) {
    bool showGui = true;

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("H2sm444.wf");
    auto ac = electronicWaveFunctionProblem.getAtomCollection();
    std::cout << ac << std::endl;


    CollectionParser collectionParser;
    auto ec = collectionParser.electronCollectionFromJson("H2sm444_TS_NRopt.json");
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


    Eigen::MatrixXd hess(n, n);
    electronicWaveFunctionProblem.hessian(x, hess);
    std::cout << hess << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> selfAdjointEigenSolver(hess,Eigen::ComputeEigenvectors);

    auto eigenvalues = selfAdjointEigenSolver.eigenvalues();
    std::cout << eigenvalues << std::endl;
    auto eigenvectors = selfAdjointEigenSolver.eigenvectors();
    std::cout << eigenvectors << std::endl;
    std::cout << std::endl;

    if (showGui) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC, "C");

        // Visualization
        MoleculeWidget moleculeWidget;
        Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

        AtomCollection3D(root, ElectronicWaveFunction::getInstance().getAtomCollection());


        // Plot eigenvectors

        int evIndex = 0;
        //for (int evIndex  = 0; evIndex  < ec.numberOfParticles(); ++evIndex ) {
            for (int i = 0; i <
                            ec.numberOfParticles(); ++i) {//for (auto i : electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei()){
                QVector3D v1(x(i * 3 + 0), x(i * 3 + 1), x(i * 3 + 2));
                QVector3D v2 = v1;
                QVector3D ev(eigenvectors.col(evIndex)(i * 3 + 0),
                             eigenvectors.col(evIndex)(i * 3 + 1),
                             eigenvectors.col(evIndex)(i * 3 + 2));
                v2 += ev*0.1;
                std::vector<QVector3D> points = {v1, v2};
                Polyline pl(root, QColor(Qt::black), points, 0.01, true);
            }
        //}

        // Plot the final point
        ElectronCollection3D(root, ElectronCollection(x, ec.spinTypesAsEigenVector()), false);

        return app.exec();
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
