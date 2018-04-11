//
// Created by Michael Heuer on 10.04.18.
//
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <QApplication>
#include <QtWidgets>

#include "ElectronicWaveFunctionProblem.h"
#include "CollectionParser.h"
#include "MoleculeWidget.h"
#include "ElectronsVector3D.h"
#include "AtomsVector3D.h"
#include "LocalNewtonSearch.h"
#include "Visualization.h"
#include "AmolqcFileImport/RefFileImporter.h"
#include "PositionsVectorTransformer.h"



Eigen::VectorXd permutePositionsCyclic(const Eigen::VectorXd &x, std::vector<unsigned> order) {
    assert(order.size() > 0);
    assert(x.size()%3 == 0);
    assert(order.size() <= x.size()/3);

    auto xnew = x;

    for (int i = 0; i < order.size()-1; ++i)
        xnew.segment(order[i]*3,3) = xnew.segment(order[i+1]*3,3);
    xnew.segment( order[order.size()-1]*3,3) = x.segment(order[0]*3,3);

    return xnew;
}


nlohmann::json eigenvectorsToJsonArray(Eigen::MatrixXd eigenvectors){

    auto jarr = nlohmann::json::array();

    for (int k = 0; k < eigenvectors.cols(); ++k) {
        auto eigenvector = nlohmann::json::array();
        auto eigenvectorFromMatrixColumn = eigenvectors.col(k);
        for (int i = 0; i < eigenvectors.rows()/3; ++i) {
            eigenvector.emplace_back(eigenvectorFromMatrixColumn[i*3+0]);
            eigenvector.emplace_back(eigenvectorFromMatrixColumn[i*3+1]);
            eigenvector.emplace_back(eigenvectorFromMatrixColumn[i*3+2]);
        }
        jarr.emplace_back(eigenvector);
    }

    return jarr;
}

nlohmann::json eigenvaluesToJsonArray(Eigen::VectorXd eigenvalues){

    auto eigarr = nlohmann::json::array();

    for (int k = 0; k < eigenvalues.size(); ++k) {
        eigarr.emplace_back(eigenvalues[k]);
    }
    return eigarr;
}


nlohmann::json resultToJSON(Eigen::MatrixXd eigenVectors) {
}

int main(int argc, char *argv[]) {

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("H6TS_CAS23_Ic444.wf", true, true);
    auto av = electronicWaveFunctionProblem.getAtomsVector();
    std::cout << av << std::endl;

    CollectionParser collectionParser;
    //RefFileImporter refFileImporter("H6TS_CAS23_Ic444.ref");
    //auto ev = refFileImporter.getMaximaStructure(1,1);
    //nlohmann::json json = collectionParser.atomsAndElectronsVectorToJson(av,ev);
    //collectionParser.writeJSON(json,"H6TS_CAS23_Ic444_GlobMax.json");

    auto ec = collectionParser.electronsVectorFromJson(collectionParser.readJSON("H6TS_CAS23_Ic444_GlobMax.json"));
    auto x1 = ec.positionsVector().positionsAsEigenVector();

    // cyclic permutation
    //auto x2 = permutePositionsCyclic(x1, {0, 3, 1, 4, 2, 5});

    // pairwise swap
    auto x2 = ec.positionsVector().positionsAsEigenVector();
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
    auto gradev = ElectronsVector(grad, ec.spinTypesVector().spinTypesAsEigenVector());
    std::cout << "Gradient:\n" << gradev << std::endl;

    Eigen::MatrixXd hess(n, n);
    electronicWaveFunctionProblem.hessian(guess, hess);

    Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> eigenSolver(hess, Eigen::ComputeEigenvectors);
    std::cout << eigenSolver.eigenvalues() << std::endl;
    std::cout << std::endl;

    // save results
    ElectronsVector ev(guess,ec.spinTypesVector().spinTypesAsEigenVector());
    nlohmann::json json = collectionParser.atomsAndElectronsVectorToJson(av,ev);
    json["Value"] = electronicWaveFunctionProblem.value(guess);
    json["Gradient"] = collectionParser.electronsVectorToJson(gradev)["ElectronsVector"];
    json["ElectronsVector"]["AtNuclei"]= electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei();
    json["ElectronsVector"]["NotAtNuclei"]= electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei();
    json["Eigenvectors"] = eigenvectorsToJsonArray(eigenSolver.eigenvectors());
    json["Eigenvalues"] = eigenvaluesToJsonArray(eigenSolver.eigenvalues());
    collectionParser.writeJSON(json,"H6TS_2pairs_swap.json");


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
            ElectronsVector ecnew(guess, ec.spinTypesVector().spinTypesAsEigenVector());
            ElectronsVector3D(root, av, ecnew, true);

            //Qt3DRender::QRenderCapture capture(root);
            //auto image = capture.requestCapture()->image();
            //image.save("image.png");

            app.exec();
        }
    }
    return 0;
}
