//
// Created by Michael Heuer on 07.11.17.
//

#include <QApplication>

#include "CollectionParser.h"
#include "ElectronicWaveFunctionProblem.h"
#include "solver/bfgsnssolver.h"
#include "solver/bfgssolver.h"
#include "solver/timeintegrationsolver.h"
#include "solver/gradientdescentumrigarlimitedsteplength.h"
#include "solver/gradientdescentsolver.h"
#include "solver/gradientdescentsimplesolver.h"
#include "solver/timeintegrationumrigarsolver.h"
#include "solver/bfgsumrigarsolver.h"

#include "AtomCollection3D.h"
#include "ElectronCollection3D.h"
#include "ParticleCollectionPath3D.h"
#include "MoleculeWidget.h"

bool handleCommandlineArguments(int argc, char **argv,
                                std::string &wavefunctionFilename,
                                std::string &electronCollectionFilename,
                                bool &showGui) {
    if (argc < 3) {
        std::cout << "Usage: \n"
                  << "Argument 1: wavefunction filename (.wf)\n"
                  << "Argument 2: electron collection filename (.json)\n"
                  << "Argument 3 (Optional): display the gui (gui)" << std::endl;
        std::cout << "Ethylene-em-5.wf LD_Ethlyen_Start.json gui" << std::endl;
        return false;
    } else if (argc >= 3) {
        wavefunctionFilename = argv[1];
        electronCollectionFilename = argv[2];
        if (argc > 3) showGui = (std::string(argv[3]) == "gui");
        return true;
    }
}

int main(int argc, char *argv[]) {
    std::string wavefunctionFilename = "H2sm444.wf"; // overwrite command line
    std::string electronCollectionFilename = "H2sm444_TS_ev.json"; // overwrite command line
    bool showGui = true;

    if( wavefunctionFilename.empty() && electronCollectionFilename.empty()) {
        bool inputArgumentsFoundQ =
                handleCommandlineArguments(argc, argv, wavefunctionFilename, electronCollectionFilename, showGui);
        if(!inputArgumentsFoundQ) return 0;
    }

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wavefunctionFilename);
    CollectionParser collectionParser;
    auto ac = electronicWaveFunctionProblem.getAtomCollection();
    auto ec = collectionParser.electronCollectionFromJson(electronCollectionFilename);
    std::cout << ac << std::endl;
    std::cout << ec << std::endl;

    Eigen::VectorXd x(ec.positionCollection().positionsAsEigenVector());
    std::cout << x.transpose() << std::endl;
    Eigen::VectorXd grad(ec.numberOfEntities());
    electronicWaveFunctionProblem.putElectronsIntoNuclei(x,grad);


    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();

    cppoptlib::TimeIntegrationSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    crit.gradNorm = 1e-6;
    crit.iterations = 1000;
    solver.setStopCriteria(crit);

    solver.minimize(electronicWaveFunctionProblem, x);
    std::cout << ElectronCollection(x,ec.spinTypeCollection().spinTypesAsEigenVector())<<std::endl;



    if(showGui) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC,"C");

        // Prepare the optimization path for visualization
        auto optimizationPath = electronicWaveFunctionProblem.getOptimizationPath();
        ElectronCollections shortenedPath(optimizationPath[0]);
        unsigned long nwanted = 300;
        auto skip = 1 + (optimizationPath.numberOfEntities() / nwanted);
        std::cout << "displaying structures with a spacing of " << skip << "." << std::endl;
        for (unsigned long i = 0; i < optimizationPath.numberOfEntities(); i = i + skip) {
            shortenedPath.append(optimizationPath[i]);
        }
        auto ecEnd = ElectronCollection(x, optimizationPath.spinTypeCollection().spinTypesAsEigenVector());
        shortenedPath.append(ecEnd);

        // Visualization
        MoleculeWidget moleculeWidget;
        Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

        AtomCollection3D(root, ElectronicWaveFunction::getInstance().getAtomCollection());

        // Plot the starting point
        ElectronCollection3D(root, ElectronCollection(x, optimizationPath.spinTypeCollection().spinTypesAsEigenVector()), false);

        // Plot the optimization path
        ParticleCollectionPath3D(root, shortenedPath);

        return app.exec();
    }
    return 0;
};
