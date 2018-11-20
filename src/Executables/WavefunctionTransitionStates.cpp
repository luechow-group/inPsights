//
// Created by Michael Heuer on 21.02.18.
//

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <QApplication>

#include "ElectronicWaveFunctionProblem.h"
#include "Serialization.h"
#include "MoleculeWidget.h"
#include "ElectronsVector3D.h"
#include "AtomsVector3D.h"
#include "LocalNewtonSearch.h"
#include "Visualization.h"
#include "../AmolqcInterface/AmolqcFileImport/include/RefFileImporter.h"
#include "PositionsVectorTransformer.h"
#include <NaturalConstants.h>

int main(int argc, char *argv[]) {

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("CP+.wf", true, true);
    auto av = electronicWaveFunctionProblem.getAtomsVector();
    std::cout << av << std::endl;


    RefFileImporter refFileImporter("CP+.ref");
    //auto ev = refFileImporter.getMaximaStructure(1,1);
    //nlohmann::json json = CollectionParser::atomsAndElectronsVectorToJson(av,ev);
    //CollectionParser::writeJSON(json,"CP+_Max1A.json")
    auto ev = YAML::LoadFile("CP+_Max1A.json")["Electrons"].as<ElectronsVector>();

    //auto ev = CollectionParser::electronsVectorFromJson(CollectionParser::readJSON("CP+_Max2A.json"));
    PositionsVector pv = ev.positionsVector();

    //auto ev =

    // MaxA1
    /*
    PositionsVector pcsel1;
    pcsel1.append(pv[ 8]);
    pcsel1.append(pv[17]);
    pcsel1.append(pv[18]);
    pcsel1.append(pv[13]);
    PositionsVector pcsel2;
    pcsel2.append(pv[ 4]);
    pcsel2.append(pv[15]);
    pcsel2.append(pv[16]);
    pcsel2.append(pv[14]);

    PositionsVectorTransformer::rotateAroundAxis(pcsel1, 60.*ConversionFactors::deg2rad,
                                                 av[0].position(), av[3].position());
    PositionsVectorTransformer::rotateAroundAxis(pcsel2, 60.*ConversionFactors::deg2rad,
                                                 av[1].position(), av[4].position());
    pv( 8) = pcsel1[0];
    pv(17) = pcsel1[1];
    pv(18) = pcsel1[2];
    pv(13) = pcsel1[3];
    pv( 4) = pcsel2[0];
    pv(15) = pcsel2[1];
    pv(16) = pcsel2[2];
    pv(14) = pcsel2[3];*/

    PositionsVector pcsel1;
    pcsel1.append(pv[15]);
    pcsel1.append(pv[17]);
    pcsel1.append(pv[18]);
    pcsel1.append(pv[ 6]);
    PositionsVector pcsel2;
    pcsel2.append(pv[ 4]);
    pcsel2.append(pv[ 8]);
    pcsel2.append(pv[ 9]);
    pcsel2.append(pv[14]);

    PositionsVectorTransformer::rotateAroundAxis(pcsel1, 60.*ConversionFactors::deg2rad,
                                                 av[0].position(), av[3].position());
    PositionsVectorTransformer::rotateAroundAxis(pcsel2,-60.*ConversionFactors::deg2rad,
                                                 av[1].position(), av[4].position());

    pv.entity(15).dataRef() = pcsel1[0];
    pv.entity(17).dataRef() = pcsel1[1];
    pv.entity(18).dataRef() = pcsel1[2];
    pv.entity( 6).dataRef() = pcsel1[3];
    pv.entity( 4).dataRef() = pcsel2[0];
    pv.entity( 8).dataRef() = pcsel2[1];
    pv.entity( 9).dataRef() = pcsel2[2];
    pv.entity(14).dataRef() = pcsel2[3];

    Eigen::VectorXd guess = pv.asEigenVector();

    // optimize the guess
    LocalNewtonSearch localNewton;
    localNewton.search(electronicWaveFunctionProblem,guess);

    auto n = ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3;
    Eigen::VectorXd grad(n);
    electronicWaveFunctionProblem.putElectronsIntoNuclei(guess,grad);



    std::cout << "gradient:" << std::endl;
    electronicWaveFunctionProblem.gradient(guess,grad);
    ElectronsVector gradev(PositionsVector(grad),ev.typesVector());
    std::cout << gradev << std::endl;

    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;
    for (auto & it : electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei()) std::cout << it << " ";
    std::cout << std::endl;


    Eigen::MatrixXd hess(n, n);
    electronicWaveFunctionProblem.hessian(guess, hess);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(hess, Eigen::ComputeEigenvectors);
    //std::cout << hess << std::endl;
    std::cout << eigenSolver.eigenvalues().transpose() << std::endl;
    std::cout << std::endl;


    // save results
    //nlohmann::json json = CollectionParser::atomsAndElectronsVectorToJson(av, ElectronsVector(PositionsVector(guess),ev.typesVector()));

    YAML::Node node;




    node["Value"] = electronicWaveFunctionProblem.value(guess);
    node["Gradient"] = gradev;
    node["ElectronsVector"]["AtNuclei"]= electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei();
    node["ElectronsVector"]["NotAtNuclei"]= electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei();
    //node["HessianDiagonalization"] = CollectionParser::selfAdjointEigenSolverResultsToJsonArray(eigenSolver);//TODO

    std::string filename = "CP+_Max2AB_TSguess.json";
    std::ofstream file(filename);
    file << Serialization::jsonStringFrom(node);


    if (true) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC, "C");

        // Visualization
        auto moleculeWidget = new MoleculeWidget();
        auto root = moleculeWidget->getRoot();

        // Plot atoms
        AtomsVector3D(root, av);

        // Plot eigenvectors
        int evIndex = 1;
        Visualization::drawEigenVector(root,eigenSolver.eigenvectors(),guess,evIndex);

        // Plot electrons with connections
        ElectronsVector ecnew(PositionsVector(guess),ev.typesVector());
        ElectronsVector3D(root, av, ecnew, true);

        return app.exec();
    }
}
