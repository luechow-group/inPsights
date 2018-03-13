#include <iostream>
#include <Eigen/Core>

#include <Qt3DCore>
#include <Qt3DRender>
#include <Qt3DExtras>

#include <QtWidgets/QApplication>

#include "MoleculeWidget.h"
#include "AtomsVector3D.h"
#include "ElectronsVector3D.h"

#include "ParticlesVectorPath3D.h"

#include "WfFileImporter.h"

#include "ElectronicWaveFunctionProblem.h"
#include <Eigen/Eigenvalues>

#include "solver/bfgssolver.h"
#include "CollectionParser.h"

#include "PotentialProblem.h"
#include "LagrangeProblem.h"

int main(int argc, char *argv[]) {
    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("Ethylene-em-5.wf");
    AtomsVector nuclei = electronicWaveFunctionProblem.getAtomsVector();

    PotentialProblem potentialProblem(nuclei);

    CollectionParser collectionParser;
    ElectronsVector electrons = collectionParser.electronsVectorFromJson("LD_Ethlyen_Start.json");

    double energy = -75;

    LagrangeProblem<PotentialProblem,ElectronicWaveFunctionProblem>
            lagrangeProblem(potentialProblem,electronicWaveFunctionProblem,energy);

    double lambdaInit = 1;

    Eigen::VectorXd x(electrons.positionsVector().positionsAsEigenVector());

    Eigen::VectorXd y = x;
    y.conservativeResize(y.size()+1,Eigen::NoChange);
    y[y.size() - 1] = lambdaInit;

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();

    cppoptlib::BfgsSolver<LagrangeProblem<PotentialProblem,ElectronicWaveFunctionProblem>> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);

    //solver.minimize(lagrangeProblem, y);

    cppoptlib::BfgsSolver<ElectronicWaveFunctionProblem> solver2;
    solver2.setDebug(cppoptlib::DebugLevel::High);
    solver2.setStopCriteria(crit);

    solver2.minimize(electronicWaveFunctionProblem,x);
    std::cout << ElectronsVector(x,electrons.spinTypesVector().spinTypesAsEigenVector())<<std::endl;

    bool showGui = true;

    if(showGui) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC,"C");

        // Prepare the optimization path for visualization
        auto optimizationPath = electronicWaveFunctionProblem.getOptimizationPath();
        ElectronsVectorCollection shortenedPath(optimizationPath[0]);
        unsigned long nwanted = 300;
        auto skip = 1 + (optimizationPath.numberOfEntities() / nwanted);
        std::cout << "displaying structures with a spacing of " << skip << "." << std::endl;
        for (unsigned long i = 0; i < optimizationPath.numberOfEntities(); i = i + skip) {
            shortenedPath.append(optimizationPath[i]);
        }
        auto ecEnd = ElectronsVector(x, optimizationPath.spinTypesVector().spinTypesAsEigenVector());
        shortenedPath.append(ecEnd);

        // Visualization
        MoleculeWidget moleculeWidget;
        Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

        AtomsVector3D(root, ElectronicWaveFunction::getInstance().getAtomsVector());

        // Plot the starting point
        ElectronsVector3D(root, ElectronsVector(x, optimizationPath.spinTypesVector().spinTypesAsEigenVector()), false);

        // Plot the optimization path
        ParticlesVectorPath3D(root, shortenedPath);

        return app.exec();
    }

    return 0;
}
