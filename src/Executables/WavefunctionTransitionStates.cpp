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
#include "PositionsVectorTransformer.h"
#include <NaturalConstants.h>

int main(int argc, char *argv[]) {

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("CP+.wf", true, true);
    auto av = electronicWaveFunctionProblem.getAtomsVector();
    std::cout << av << std::endl;


    RefFileImporter refFileImporter("CP+.ref");
    auto ev = refFileImporter.getMaximaStructure(1,1);
    nlohmann::json json = CollectionParser::atomsAndElectronsVectorToJson(av,ev);
    CollectionParser::writeJSON(json,"CP+_Guess.json");

    auto ec = CollectionParser::electronsVectorFromJson(CollectionParser::readJSON("CP+_Guess.json"));
    PositionsVector pc = ec.positionsVector();

    PositionsVector pcsel1;
    pcsel1.append(pc[8]);
    pcsel1.append(pc[17]);
    pcsel1.append(pc[18]);
    pcsel1.append(pc[13]);

    PositionsVector pcsel2;
    pcsel2.append(pc[4]);
    pcsel2.append(pc[15]);
    pcsel2.append(pc[16]);
    pcsel2.append(pc[14]);

    //std::cout << pcsel2 << std::endl;
    PositionsVectorTransformer::rotateAroundAxis(pcsel1, 120.*ConversionFactors::deg2rad,
                                                 av[0].position(), av[3].position());
    PositionsVectorTransformer::rotateAroundAxis(pcsel2, 120.*ConversionFactors::deg2rad,
                                                 av[1].position(), av[4].position());
    //auto x2 = x1;

    pc( 8) = pcsel1[0];
    pc(17) = pcsel1[1];
    pc(18) = pcsel1[2];
    pc(13) = pcsel1[3];
    pc( 4) = pcsel2[0];
    pc(15) = pcsel2[1];
    pc(16) = pcsel2[2];
    pc(14) = pcsel2[3];

    Eigen::VectorXd guess = pc.positionsAsEigenVector();

    // optimize the guess
    LocalNewtonSearch localNewton;
    localNewton.search(electronicWaveFunctionProblem,guess);

    auto n = ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3;
    Eigen::VectorXd grad(n);
    electronicWaveFunctionProblem.putElectronsIntoNuclei(guess,grad);
    electronicWaveFunctionProblem.gradient(guess,grad);


    std::cout << "gradient:" << std::endl;
    std::cout << ElectronsVector(grad,ec.spinTypesVector().spinTypesAsEigenVector()) << std::endl;

    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;
    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;


    Eigen::MatrixXd hess(n, n);
    electronicWaveFunctionProblem.hessian(guess, hess);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(hess, Eigen::ComputeEigenvectors);
    auto eigenvectors = eigenSolver.eigenvectors();
    std::cout << hess << std::endl;
    std::cout << eigenSolver.eigenvalues() << std::endl;
    std::cout << std::endl;


    //// save results
    //ElectronsVector ev(guess,ec.spinTypesVector().spinTypesAsEigenVector());
    //nlohmann::json json = CollectionParser::atomsAndElectronsVectorToJson(av,ev);
    //json["Value"] = electronicWaveFunctionProblem.value(guess);
    //json["Gradient"] = CollectionParser::electronsVectorToJson(gradev)["ElectronsVector"];
    //json["ElectronsVector"]["AtNuclei"]= electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei();
    //json["ElectronsVector"]["NotAtNuclei"]= electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei();
    //json["HessianDiagonalization"] = CollectionParser::selfAdjointEigenSolverResultsToJsonArray(eigenSolver);
    //CollectionParser::writeJSON(json,"H6TS_2pairs_swap.json");



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
        Visualization::drawEigenVector(root,eigenvectors,guess,evIndex);

        // Plot electrons with connections
        ElectronsVector ecnew(guess,ec.spinTypesVector().spinTypesAsEigenVector());
        ElectronsVector3D(root, av, ecnew, true);

        return app.exec();
    }
}
