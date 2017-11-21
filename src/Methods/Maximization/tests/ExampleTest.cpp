//
// Created by dahl on 21.11.17.
//

#include <gmock/gmock.h>
#include <Eigen/Core>
#include <iostream>
#include "ElectronicWaveFunctionProblem.h"
#include "solver/gradientdescentumrigarlimitedsteplength.h"

using namespace testing;

class AGradientDescentUmrigarLimitedStepLengthSolverTest : public Test {};

TEST_F(AGradientDescentUmrigarLimitedStepLengthSolverTest , Ethane) {

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
    cppoptlib::GradientDescentUmrigarLimitedSteplength<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);
    solver.setMaxStepLength(1e-1);
    solver.setSteepestDescentRate(1.0);
    solver.setDistanceCriteriaUmrigar(0.5);
    solver.setThreshholdUmrigar(1e-5);
    solver.minimize(electronicWaveFunctionProblem, xA);

    ASSERT_TRUE(false);

}