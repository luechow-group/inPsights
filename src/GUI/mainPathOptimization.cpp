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
#include "solver/timeintegrationumrigarsolver.h"
#include "solver/bfgsumrigarsolver.h"

#include "AtomCollection3D.h"
#include "ElectronCollection3D.h"
#include "ParticleCollectionPath3D.h"
#include "MoleculeWidget.h"

int main(int argc, char *argv[]) {

    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("Ethylene-em-5.wf");

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();

    /*
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

    /*
    Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);

    // For diborane: Diborane.wf; Pfad 1, permutiert
    xA << \
    -1.8108290, -0.0103127, -0.0733318,\
    -1.0669983, -0.4967279, -1.4553791,\
    1.6199701, -1.5930971, -0.2339874,\
   -2.7480445, 2.0910952, 0.0152669,\
    3.4598213, 0.7427135, 0.9798742,\
    -0.1999113, 0.0572744, 2.2458029,\
    1.5627554, -0.1209528, 0.1947172,\
    -3.4694941, -2.6586327, 0.5373582,\
    2.4173843, 0.0920298, 0.1590485,\
    2.7100060, -1.9063833, -0.1915220,\
    -2.4599213, 0.4676680, 0.5608425,\
    0.9639635, 0.9456205, 0.7843612,\
    -0.3894725, -0.2019606, -1.2126405,\
    1.6038582, 0.5784657, -0.5778354,\
    -1.6236944, -1.0465345, 0.6120425,\
    -1.6807995, -0.2473476, 0.0670798;
    */

    /*
    Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);

    // For fluoride: Fluoride.wf
    xA << \
    0.714583,  2.171709,  2.377429,\
   -0.805267,  0.373607,  0.961730,\
   -0.013201, -0.133104,  1.591434,\
    0.305421,  0.141948, -1.348573,\
    0.733193,  0.981820, -1.175646,\
   -1.631408,  0.621465, -2.145863,\
   -0.549092, -2.641827,  3.075085,\
    1.668276,  1.311450, -0.564745,\
    1.554871, -1.617958, -2.978076,\
    0.204682, -0.074675,  1.290531,\
   -0.880830, -0.567550,  0.091367,\
    0.483278, -2.256104,  1.174406,\
    1.888764, -0.491579, -1.046459,\
   -0.154281,  1.014234, -2.217571,\
   -0.924312,  0.945934, -0.019794,\
   -1.987497,  0.072370,  1.736939,\
   -0.544636, -2.204059, -3.499582,\
    0.005195,  0.207915, -1.906905;
    */

    /*
    Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);

    // For H2O2: h2o2_rb3_tzp_em.wf
    xA << \
    1.955104, 0.926238, -0.398171,\
    -0.955123, 0.326257,0.097756,\
    1.740730, -0.563830, 0.983537,\
    2.007865, 0.350172, 0.564486,\
    0.116030, 0.305585, -0.810662,\
    -2.428552, -0.199389, -0.407413,\
    -1.900057, -0.994862, 0.990089,\
    1.290431, 0.709741, -0.676654,\
    -1.917662, -0.109087, -0.231133,\
    1.895104, 0.566238, -0.148171,\
    -1.155123, -0.006257, -0.357756,\
    -2.540503, 1.553981, 1.103858,\
    -2.541622, -1.716824, 0.821543,\
    1.493530, -0.204730, -0.305049,\
    2.200479, 0.350482, 0.753813,\
    -1.197279, -1.009831, 0.031789,\
    -2.288358, 0.393611, -1.051923,\
    0.574047, 1.639305, -0.509279;
    */


    Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);

    // For Ethylen: Ethylene-em-5.wf
    xA << \
    0.714583,  2.171709,  2.377429,\
   -0.805267,  0.373607,  0.961730,\
   -0.013201, -0.133104,  1.591434,\
    0.305421,  0.141948, -1.348573,\
    0.733193,  0.981820, -1.175646,\
   -1.631408,  0.621465, -2.145863,\
   -0.549092, -2.641827,  3.075085,\
    1.668276,  1.311450, -0.564745,\
    1.554871, -1.617958, -2.978076,\
    0.204682, -0.074675,  1.290531,\
   -0.880830, -0.567550,  0.091367,\
    0.483278, -2.256104,  1.174406,\
    1.888764, -0.491579, -1.046459,\
   -0.154281,  1.014234, -2.217571,\
   -0.924312,  0.945934, -0.019794,\
   -1.987497,  0.072370,  1.736939;
   //-0.544636, -2.204059, -3.499582,\
   // 0.005195,  0.207915, -1.906905;


    cppoptlib::BfgsUmrigarSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    crit.gradNorm = 1e-5;
    crit.iterations = 100000;
    solver.setStopCriteria(crit);
    //solver.setMaxStepLength(1.0);
    //solver.setSteepestDescentRate(0.1);
    //solver.setDistanceCriteriaUmrigar(0.1);

    solver.minimize(electronicWaveFunctionProblem, xA);

    std::cout<<xA<<std::endl;

    auto optimizationPath = electronicWaveFunctionProblem.getOptimizationPath();
    ElectronCollections shortenedPath(ElectronCollection(optimizationPath.front(),
                                                         optimizationPath.getSpinTypeCollection()));

    unsigned long nwanted = 300;
    auto skip = 1+(optimizationPath.length()/nwanted);
    std::cout << "displaying structures with a spacing of " << skip << "." << std::endl;
    for (unsigned long i = 0; i < optimizationPath.length(); i=i+skip) {
        shortenedPath.append(ElectronCollection(optimizationPath.getElectronCollection(i),
                                                optimizationPath.getSpinTypeCollection()));
    }

    shortenedPath.append(ElectronCollection(xA,optimizationPath.getSpinTypeCollection().spinTypesAsEigenVector()));

    //visualization
    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    AtomCollection3D(root,ElectronicWaveFunction::getInstance().getAtomCollection());

    // Plot the starting point
    ElectronCollection3D(root, ElectronCollection(ParticleCollection(xA),
                                                  optimizationPath.getSpinTypeCollection()), false);

    // Plot the optimization path
    ParticleCollectionPath3D(root, shortenedPath);

    return app.exec();
};

