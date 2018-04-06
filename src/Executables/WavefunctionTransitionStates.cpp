//
// Created by Michael Heuer on 21.02.18.
//

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <QApplication>

#include "ElectronicWaveFunctionProblem.h"
#include "CollectionParser.h"
#include "MoleculeWidget.h"
#include "ElectronsVector3D.h"
#include "AtomsVector3D.h"
#include "Polyline.h"
#include "LocalNewtonSearch.h"

int main(int argc, char *argv[]) {

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("H2sm444.wf");
    auto ac = electronicWaveFunctionProblem.getAtomsVector();
    std::cout << ac << std::endl;
    
    CollectionParser collectionParser;
    auto ec = collectionParser.electronsVectorFromJson("H2sm444_TSguess.json");
    auto x = ec.positionsVector().positionsAsEigenVector();
    std::cout << ec << std::endl;

    // optimize the guess
    LocalNewtonSearch localNewton;
    localNewton.search(electronicWaveFunctionProblem,x);
    ec.positionsVector() = PositionsVector(x);


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

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(hess, Eigen::ComputeEigenvectors);
    auto eigenvectors = eigenSolver.eigenvectors();
    std::cout << hess << std::endl;
    std::cout << eigenSolver.eigenvalues() << std::endl;
    std::cout << std::endl;


    if (true) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC, "C");

        // Visualization
        MoleculeWidget moleculeWidget;
        Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

        AtomsVector3D(root, ElectronicWaveFunction::getInstance().getAtomsVector());


        // Plot eigenvectors
        int evIndex = 0;
        for (int i = 0; i < ec.numberOfEntities(); ++i) {
            QVector3D v1(x(i * 3 + 0), x(i * 3 + 1), x(i * 3 + 2));
            QVector3D v2 = v1;
            QVector3D ev(eigenvectors.col(evIndex)(i * 3 + 0),
                         eigenvectors.col(evIndex)(i * 3 + 1),
                         eigenvectors.col(evIndex)(i * 3 + 2));
            v2 += ev;
            std::vector<QVector3D> points = {v1, v2};
            Polyline pl(root, QColor(Qt::black), points, 0.01, true);
        }

        // Plot the final point
        ElectronsVector3D(root, ec, false);

        return app.exec();
    }
}
