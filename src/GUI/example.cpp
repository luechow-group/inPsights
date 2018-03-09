//
// Created by Michael Heuer on 07.11.17.
//

#include <QApplication>
#include <Qt3DRender>
#include <Qt3DExtras>
#include <QtWidgets>
#include <iostream>

#include "CollectionParser.h"
#include "OptimizationPathFileImporter.h"
#include "ElectronicWaveFunctionProblem.h"
#include "solver/bfgsnssolver.h"
#include "solver/bfgssolver.h"
#include "solver/timeintegrationsolver.h"
#include "solver/gradientdescentumrigarlimitedsteplength.h"
#include "solver/gradientdescentsolver.h"
#include "solver/gradientdescentsimplesolver.h"

#include "AtomsVector3D.h"
#include "ElectronsVector3D.h"
#include "ParticlesVectorPath3D.h"
#include "MoleculeWidget.h"

int main(int argc, char *argv[]) {

    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("Ethane-em-5.wf");

    CollectionParser collectionParser;
    auto ecA = collectionParser.electronsVectorFromJson("Ethane-glob-max.json");
    auto xA = ecA.positionsAsEigenVector();

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();

    //OptimizationPathFileImporter optimizationPathFileImporter("Diborane-Paths.300",1); // Aufpassen ob richtige Multiplizität
    //std::string wfFilename = "Diborane.wf";
    //ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wfFilename);
    //auto numberOfPaths = optimizationPathFileImporter.getNumberOfPaths();
    //auto psiSquareDistributedParticlesVector = optimizationPathFileImporter.getPath(6).front();
    //VectorXd xA = psiSquareDistributedParticlesVector.positionsAsEigenVector();

    cppoptlib::TimeIntegrationSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    crit.gradNorm = 1e-5;
    crit.iterations = 1000;
    solver.setStopCriteria(crit);
    //solver.setMaxStepLength(0.2);
    //solver.setSteepestDescentRate(0.2);
    //solver.setDistanceCriteriaUmrigar(0.5);

    solver.minimize(electronicWaveFunctionProblem, xA);

    // export result to json
    auto ecAresult = ElectronsVector(xA,ecA.spinTypesAsEigenVector());
    auto jsonObject = collectionParser.electronsVectorToJson(ecAresult);
    jsonObject["comment"]= "optimized with FIRE";

    nlohmann::json solverSettings;
    solverSettings["iterations"] = crit.iterations;
    solverSettings["gradNorm"] = crit.gradNorm;
    jsonObject["settings"] = solverSettings;

    // write to file
    collectionParser.writeJSON(jsonObject,"FIRE-Ehane-opt.json");

    auto optimizationPath = electronicWaveFunctionProblem.getOptimizationPath();
    ElectronsVectorCollection shortenedPath(ElectronsVector(optimizationPath.front(),
                                                         optimizationPath.getSpinTypesVector()));

    unsigned long nwanted = 300;
    auto skip = 1+(optimizationPath.length()/nwanted);
    std::cout << "displaying structures with a spacing of " << skip << "." << std::endl;
    for (unsigned long i = 0; i < optimizationPath.length(); i=i+skip) {
        shortenedPath.append(ElectronsVector(optimizationPath.getElectronsVector(i),
                                                optimizationPath.getSpinTypesVector()));
    }

    shortenedPath.append(ElectronsVector(xA,optimizationPath.getSpinTypesVector().spinTypesAsEigenVector()));

    //visualization
    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    AtomsVector3D(root,ElectronicWaveFunction::getInstance().getAtomsVector());

    // Plot the starting point
    ElectronsVector3D(root, ElectronsVector(ParticlesVector(xA),
                                                  optimizationPath.getSpinTypesVector()), true);

    // Plot the optimization path
    ParticlesVectorPath3D(root, shortenedPath);

    return app.exec();
};