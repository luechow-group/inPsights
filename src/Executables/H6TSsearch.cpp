//
// Created by Michael Heuer on 10.04.18.
//

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <QApplication>
#include <QtWidgets>

#include "ElectronicWaveFunctionProblem.h"
#include "AmolqcFileImport/RefFileImporter.h"
#include "Serialization.h"
#include "PositionsVectorTransformer.h"
#include "MoleculeWidget.h"
#include "ElectronsVector3D.h"
#include "AtomsVector3D.h"
#include "Visualization.h"
#include "LocalNewtonSearch.h"




int main(int argc, char *argv[]) {

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("H6TS_CAS23_Ic444.wf", true, true);
    auto av = electronicWaveFunctionProblem.getAtomsVector();
    std::cout << av << std::endl;

    //RefFileImporter refFileImporter("H6TS_CAS23_Ic444.ref");
    //auto ev = refFileImporter.getMaximaStructure(1,1);
    //nlohmann::json json = CollectionParser::atomsAndElectronsVectorToJson(av,ev);
    //CollectionParser::writeJSON(json,"H6TS_CAS23_Ic444_GlobMax.json");

    auto ec = YAML::LoadFile("H6TS_CAS23_Ic444_GlobMax.json")["Electrons"].as<ElectronsVector>();
    auto x1 = ec.positionsVector().asEigenVector();

    // cyclic permutation
    //auto x2 = PositionsVectorTransfomer::permutePositionsCyclic(x1, {0, 3, 1, 4, 2, 5});

    // pairwise swap
    auto x2 = ec.positionsVector().asEigenVector();
    x2.segment(0*3,3) = x1.segment(3*3,3);
    x2.segment(3*3,3) = x1.segment(0*3,3);
    x2.segment(1*3,3) = x1.segment(4*3,3);
    x2.segment(4*3,3) = x1.segment(1*3,3);
    x2.segment(2*3,3) = x1.segment(5*3,3);
    x2.segment(5*3,3) = x1.segment(2*3,3);

    Eigen::VectorXd guess = 0.5 * (x2 - x1) + x1;
    // z-deviation
    guess(0*3+2) += 0.2;
    guess(3*3+2) -= 0.2;
    guess(1*3+2) -= 0.2;
    guess(4*3+2) += 0.2;
    guess(2*3+2) += 0.2;
    guess(5*3+2) -= 0.2;

    //auto noise = Eigen::VectorXd::Random(guess.size());
    //guess += 0.05*noise;



    // optimize the guess
    LocalNewtonSearch localNewton;
    localNewton.search(electronicWaveFunctionProblem, guess);

    auto n = ElectronicWaveFunction::getInstance().getNumberOfElectrons() * 3;
    Eigen::VectorXd grad(n);
    electronicWaveFunctionProblem.putElectronsIntoNuclei(guess, grad);
    electronicWaveFunctionProblem.gradient(guess, grad);
    auto gradev = ElectronsVector(PositionsVector(grad),ec.typesVector());
    std::cout << "Gradient:\n" << gradev << std::endl;

    Eigen::MatrixXd hess(n, n);
    electronicWaveFunctionProblem.hessian(guess, hess);

    Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> eigenSolver(hess, Eigen::ComputeEigenvectors);
    std::cout << eigenSolver.eigenvalues() << std::endl;
    std::cout << std::endl;

    // save results
    //ElectronsVector ev(guess,ec.typesVector().asEigenVector());
    //nlohmann::json json = CollectionParser::atomsAndElectronsVectorToJson(av,ev);
    //json["Value"] = electronicWaveFunctionProblem.value(guess);
    //json["Gradient"] = CollectionParser::electronsVectorToJson(gradev)["ElectronsVector"];
    //json["ElectronsVector"]["AtNuclei"]= electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei();
    //json["ElectronsVector"]["NotAtNuclei"]= electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei();
    ////json["HessianDiagonalization"] = CollectionParser::selfAdjointEigenSolverResultsToJsonArray(eigenSolver);
    //CollectionParser::writeJSON(json,"H6TS_2pairs_swap.json");


    for (int i = 0; i < 17; ++i) {

        if (true) {
            QApplication app(argc, argv);
            setlocale(LC_NUMERIC, "C");

            // Visualization
            MoleculeWidget moleculeWidget;
            Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

            // Plot atoms
            AtomsVector3D(root, av);

            // Plot eigenvectors
            int evIndex = i;
            Visualization::drawEigenVector(root, eigenSolver.eigenvectors(), guess, evIndex);

            // Plot electrons with connections
            ElectronsVector ecnew(PositionsVector(guess), ec.typesVector());
            ElectronsVector3D(root, av, ecnew, true);

            //Qt3DRender::QRenderCapture capture(root);
            //auto image = capture.requestCapture()->image();
            //image.save("image.png");

            app.exec();
        }
    }
    return 0;
}
