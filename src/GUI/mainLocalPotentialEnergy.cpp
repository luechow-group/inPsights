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
#include "solver/gradientdescentsolver.h"
#include "solver/gradientdescentsimplesolver.h"
#include "CollectionParser.h"

#include "PotentialProblem.h"
#include "LagrangeProblem.h"
#include "GradientSqMagnitudeProblem.h"

int visualizeOptPath(const QApplication &app,
                     const ElectronsVectorCollection &optimizationPath,
                     const unsigned long &nwanted);

using namespace Eigen;

int main(int argc, char *argv[]) {
    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("H2.wf");
    AtomsVector nuclei = electronicWaveFunctionProblem.getAtomsVector();

    PotentialProblem potentialProblem(nuclei);

    CollectionParser collectionParser;

    //ElectronsVector electrons = collectionParser.electronsVectorFromJson("LD_Diboran_Start.json");

    ElectronsVector electrons;
    electrons.append(Electron(Vector3d(0,0,0.1),Spin::SpinType::alpha));
    electrons.append(Electron(Vector3d(0,0,-0.1),Spin::SpinType::beta));

    double energy = -1.7;

    LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>
            lagrangeProblem(electronicWaveFunctionProblem,potentialProblem,energy);

    GradientSqMagnitudeProblem<LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>> mainprob(lagrangeProblem);

    double lambdaInit = -1;

    Eigen::VectorXd x(electrons.positionsVector().positionsAsEigenVector());

    Eigen::VectorXd y = x;
    y.conservativeResize(y.size()+1,Eigen::NoChange);
    y[y.size() - 1] = lambdaInit;

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
 
    cppoptlib::GradientDescentSolver<GradientSqMagnitudeProblem<LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>>> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);

    solver.minimize(mainprob, y);

    bool showGui = true;

    if(showGui) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC,"C");

        // Prepare the optimization path for visualization
        auto optimizationPath = mainprob.getProblem().getProblem().getOptimizationPath();

        //prepend starting position
        optimizationPath.prepend(electrons);

        //append end position
        optimizationPath.append(ElectronsVector(PositionsVector(y.head(y.size() - 1)),
                                                optimizationPath.spinTypesVector()));

        return visualizeOptPath(app, optimizationPath, 300);

    }

    return 0;
}

int visualizeOptPath(const QApplication &app, const ElectronsVectorCollection &optimizationPath, const unsigned long &nwanted) {
    ElectronsVectorCollection shortenedPath(optimizationPath[0]);

    double optPathLength = 0.0;
    Eigen::VectorXd pathLengthVector =
            Eigen::VectorXd::Zero(optimizationPath.numberOfEntities());

    for (unsigned long i = 1; i < optimizationPath.numberOfEntities(); i++){
        optPathLength += optimizationPath.norm(i,i - 1);
        pathLengthVector[i] = optPathLength;
    }

    double stepLength = optPathLength / nwanted;
    ElectronsVector elecVector;
    unsigned long index = 0;
    for (unsigned long i = 0; i <= nwanted; i++) {
        index = 0;
        for (unsigned long j = 1; j < optimizationPath.numberOfEntities(); j++){
            if (fabs(pathLengthVector[j] - i * stepLength) <
                    fabs(pathLengthVector[index] - i * stepLength)){
                index = j;
            }
        }
        shortenedPath.append(optimizationPath[index]);
    }

    // Visualization
    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    AtomsVector3D(root, ElectronicWaveFunction::getInstance().getAtomsVector());

    // Plot the starting point
    ElectronsVector3D(root, optimizationPath[-1], false);

    // Plot the optimization path
    ParticlesVectorPath3D(root, shortenedPath);

    return app.exec();
}
