//
// Created by Michael Heuer on 07.11.17.
//

#include <QApplication>
#include <Qt3DRender>
#include <Qt3DExtras>
#include <QtWidgets>
#include <iostream>

#include "OptimizationPathFileImporter.h"
#include "RefFileImporter.h"
#include "ElectronicWaveFunction.h"
#include "ElectronicWaveFunctionProblem.h"
#include "solver/bfgsnssolver.h"
#include "solver/timeintegrationsolver.h"
#include "solver/gradientdescentumrigarlimitedsteplength.h"
#include "solver/gradientdescentsolver.h"
#include "AtomCollection3D.h"
#include "ElectronCollection3D.h"
#include "ParticleCollectionPath3D.h"
#include "MoleculeWidget.h"

int main(int argc, char *argv[]) {

    QApplication app(argc, argv);

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("Ethane-em-5.wf");
    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();

    cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    crit.gradNorm = 1e-6;
    crit.iterations = 265;
    solver.setStopCriteria(crit);

    //RefFileImporter refFileImporter("Ethane-max.ref");
    //auto maximumOfInterest = refFileImporter.getMaximaStructure(1,1);
    //auto maximumElectronPositionVector = maximumOfInterest.positionsAsEigenVector();
    //std::cout << maximumElectronPositionVector << std::endl;

    //ElectronicWaveFunction::getInstance().setFrozenElectrons({1,2,3,4,5, 10,11,12,13,14});
    //Ethane global max
    /*Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
    xA << \
    0.000000, 0.000000, 1.443184,\
    0.000000, 0.000000,-1.443184,\
   -1.662146,-0.959641, 2.192989,\
   -1.662146, 0.959641,-2.192989,\
    0.000000,-1.919300,-2.192989,\
   -0.024099, 0.773535, 1.718336+0.3,\
    0.657845,-0.407636, 1.718335-0.3,\
    0.658424, 0.380140,-1.752435+0.3,\
   -0.034848,-0.020119,-0.660100-0.3,\
    0.000000, 0.000000, 1.443184,\
    0.000000, 0.000000,-1.443184,\
    1.662146,-0.959641, 2.192989,\
    1.662146, 0.959641,-2.192989,\
    0.000000, 1.919300, 2.192989,\
   -0.657846, 0.407635,-1.718336-0.3,\
    0.024100,-0.773537,-1.718336+0.3,\
   -0.658423,-0.380140, 1.752435-0.3,\
    0.034847, 0.020120, 0.660096+0.3;*/

    /*amolqcInput1*/
    Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
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

    solver.minimize(electronicWaveFunctionProblem, xA);
    auto optimizationPath = electronicWaveFunctionProblem.getOptimizationPath();

    
    
    //visualization
    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    AtomCollection3D(root, ElectronicWaveFunction::getInstance().getAtomCollection());

    // Plot the starting point
    ElectronCollection3D(root,ElectronCollection(ParticleCollection(xA),
                                                 optimizationPath.getSpinTypeCollection()));
    
    // Plot the optimization path
    ParticleCollectionPath3D(root, optimizationPath);

    return app.exec();
};