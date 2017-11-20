//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>

#include "solver/bfgsnssolver.h"
#include "ElectronicWaveFunction.h"
#include "ElectronicWaveFunctionProblem.h"
#include "RefFileImporter.h"

int main(int argc, char const *argv[]) {
    
    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("Ethane-em-5.wf");
    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
    crit.gradNorm = 1e-6;
    crit.iterations = 100;
    cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);

    // 
    RefFileImporter refFileImporter("Ethane-max.ref");
    auto maximumOfInterest = refFileImporter.getMaximaStructure(1,1);
    auto xA = maximumOfInterest.positionsAsEigenVector();
    std::cout << xA << std::endl;

    //ElectronicWaveFunction::getInstance().setFrozenElectrons({1,2,3,4,5, 10,11,12,13,14});

    /*amolqcInput1*/
    /*Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
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
    0.005195,  0.207915, -1.906905;*/


    solver.minimize(electronicWaveFunctionProblem, xA);

}
