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


int main(int argc, char *argv[]) {

    // handle command line arguments
    std::string wavefunctionFilename = argv[1];
    std::string electronCollectionFilename = argv[2];
    bool showGui = false;
    if (argc > 3) showGui = (std::string(argv[3]) == "gui");

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wavefunctionFilename);
    CollectionParser collectionParser;

    auto ac = electronicWaveFunctionProblem.getAtomCollection();
    auto ec = collectionParser.electronCollectionFromJson(electronCollectionFilename);
    std::cout << ac << std::endl;
    std::cout << ec << std::endl;

    Eigen::VectorXd x(ec.positionsAsEigenVector());
    std::cout << x.transpose() << std::endl;

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();

    cppoptlib::BfgsUmrigarSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    crit.gradNorm = 1e-5;
    crit.iterations = 100000;
    solver.setStopCriteria(crit);
    //solver.setMaxStepLength(1.0);
    //solver.setSteepestDescentRate(0.1);
    //solver.setDistanceCriteriaUmrigar(0.1);

    solver.minimize(electronicWaveFunctionProblem, x);
    std::cout << ElectronCollection(ParticleCollection(x),SpinTypeCollection(ec.spinTypesAsEigenVector()))<<std::endl;


    if(showGui) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC,"C");

        // Prepare the optimization path for visualization
        auto optimizationPath = electronicWaveFunctionProblem.getOptimizationPath();
        ElectronCollections shortenedPath(ElectronCollection(optimizationPath.front(),
                                                             optimizationPath.getSpinTypeCollection()));
        unsigned long nwanted = 300;
        auto skip = 1 + (optimizationPath.length() / nwanted);
        std::cout << "displaying structures with a spacing of " << skip << "." << std::endl;
        for (unsigned long i = 0; i < optimizationPath.length(); i = i + skip) {
            shortenedPath.append(ElectronCollection(optimizationPath.getElectronCollection(i),
                                                    optimizationPath.getSpinTypeCollection()));
        }

        shortenedPath.append(ElectronCollection(x, optimizationPath.getSpinTypeCollection().spinTypesAsEigenVector()));

        // Visualization
        MoleculeWidget moleculeWidget;
        Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

        AtomCollection3D(root, ElectronicWaveFunction::getInstance().getAtomCollection());

        // Plot the starting point
        ElectronCollection3D(root, ElectronCollection(ParticleCollection(x),
                                                      optimizationPath.getSpinTypeCollection()), false);

        // Plot the optimization path
        ParticleCollectionPath3D(root, shortenedPath);

        return app.exec();
    }
    return 0;
};

/* Old code from Lennard Dahl
#include "OptimizationPathFileImporter.h"
// For diborane: Diborane.wf; Pfad 6 & 7
OptimizationPathFileImporter optimizationPathFileImporter("Diborane-Paths.300",1); // Aufpassen ob richtige Multiplizität
//std::string wfFilename = "Diborane.wf";
//ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wfFilename);
//auto numberOfPaths = optimizationPathFileImporter.getNumberOfPaths();
//auto psiSquareDistributedParticleCollection = optimizationPathFileImporter.getPath(6).front();
//VectorXd xA = psiSquareDistributedParticleCollection.positionsAsEigenVector();

auto numberOfPaths = optimizationPathFileImporter.getNumberOfPaths();

auto psiSquareDistributedParticleCollection = optimizationPathFileImporter.getPath(6).front();
VectorXd xA = psiSquareDistributedParticleCollection.positionsAsEigenVector();
*/