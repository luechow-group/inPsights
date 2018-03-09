//
// Created by dahl on 21.11.17.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include <iostream>
#include <solver/gradientdescentsimplesolver.h>
#include "ElectronicWaveFunctionProblem.h"
#include "solver/gradientdescentumrigarlimitedsteplength.h"
#include "solver/bfgsnssolver.h"
#include <OptimizationPathFileImporter.h>
#include "solver/gradientdescentsolver.h"

using namespace testing;

class AGradientDescentUmrigarLimitedStepLengthSolverTest : public Test {};

    /*
    std::string wfFilename = "H2.wf";
    ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wfFilename);

    Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
    xA << \
    0.00000,  0.00000,  0.50000,\
    0.00000,  0.00000, -0.50000;*/


    //cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();

    /*
    crit.iterations = 100;
    crit.gradNorm = 1e-5;
    cppoptlib::GradientDescentUmrigarLimitedSteplength<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);
    solver.setMaxStepLength(1e-2);
    solver.setSteepestDescentRate(1.0);
    solver.setDistanceCriteriaUmrigar(0.5);
    solver.minimize(electronicWaveFunctionProblem, xA);
    */

    /*
    OptimizationPathFileImporter optimizationPathFileImporter("Diborane-Paths.300",1); // Aufpassen ob richtige Multiplizität
    std::string wfFilename = "Diborane.wf";
    ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wfFilename);

    electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei();

    auto numberOfPaths = optimizationPathFileImporter.getNumberOfPaths();

    for (unsigned long k = 1; k < numberOfPaths; ++k) {
        auto psiSquareDistributedParticlesVector = optimizationPathFileImporter.getPath(k).front();
        VectorXd x0 = psiSquareDistributedParticlesVector.positionsAsEigenVector();

        cppoptlib::GradientDescentSolver<ElectronicWaveFunctionProblem> solver;
        solver.setDebug(cppoptlib::DebugLevel::High);
        crit.gradNorm = 1e-5;
        crit.iterations = 100;
        solver.setStopCriteria(crit);
        //solver.setMaxStepLength(1e-1);
        //solver.setSteepestDescentRate(1.0);
        //solver.setDistanceCriteriaUmrigar(0.1);

        solver.minimize(electronicWaveFunctionProblem, x0);
        std::cout<<electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei().back()<<std::endl;
    }*/


    //std::cout<<x0<<std::endl;

    //ASSERT_TRUE(false);

//}

/*


    xA << \
    0.00000,  0.00000,  0.50000,\
    0.00000,  0.00000, -0.50000;

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

TEST_F(AGradientDescentUmrigarLimitedStepLengthSolverTest , 2Ethane) {

    std::string wfFilename = "Ethane-em-5.wf";
    ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wfFilename);

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


    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
    crit.iterations = 1000;
    crit.gradNorm = 1e-5;
    cppoptlib::GradientDescentSimpleSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);
    solver.minimize(electronicWaveFunctionProblem, xA);
    std::cout<<xA<<std::endl;

    ASSERT_TRUE(false);
}*/