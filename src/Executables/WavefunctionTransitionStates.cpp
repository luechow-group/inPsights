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
#include "LocalNewtonSearch.h"
#include "Visualization.h"
#include "AmolqcFileImport/RefFileImporter.h"

int main(int argc, char *argv[]) {

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("H6TS_CAS23_Ic444.wf");
    auto av = electronicWaveFunctionProblem.getAtomsVector();
    std::cout << av << std::endl;

    RefFileImporter refFileImporter("H6TS_CAS23_Ic444.ref");
    auto ev = refFileImporter.getMaximaStructure(1,1);

    CollectionParser collectionParser;

    nlohmann::json json = collectionParser.atomsAndElectronsVectorToJson(av,ev);

    collectionParser.writeJSON(json,"H6TS_CAS23_Ic444_Guess.json");

    auto stv = ev.spinTypesVector();
    auto pv = ev.positionsVector();
    PositionsVectorCollection pvc;
    pvc.append(pv);
    pvc.append(pv);
    ElectronsVectorCollection evc(pvc,stv);

    nlohmann::json json2 = collectionParser.electronsVectorCollectionToJson(evc);
    collectionParser.writeJSON(json2,"evctest.json");

    auto ec = collectionParser.electronsVectorFromJson(collectionParser.readJSON("H6TS_CAS23_Ic444_Guess.json"));
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

        // Plot atoms
        AtomsVector3D(root, av);

        // Plot eigenvectors
        int evIndex = 0;
        Visualization::drawEigenVector(root,eigenvectors,x,evIndex);

        // Plot electrons with connections
        ElectronsVector3D(root, av, ec, false);

        return app.exec();
    }
}
