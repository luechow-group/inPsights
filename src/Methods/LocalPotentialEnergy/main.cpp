//
// Created by leonard on 21.02.18.
//

#include "CollectionParser.h"
#include "ElectronicWaveFunctionProblem.h"
#include "PotentialProblem.h"
#include "LagrangeProblem.h"
#include "solver/bfgssolver.h"

int main(){
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

    /*cppoptlib::BfgsSolver<ElectronicWaveFunctionProblem> solver2;
    solver2.setDebug(cppoptlib::DebugLevel::High);
    solver2.setStopCriteria(crit);

    solver2.minimize(electronicWaveFunctionProblem,x);
    std::cout << ElectronsVector(x,electrons.spinTypesVector().spinTypesAsEigenVector())<<std::endl;*/

    return 0;
}
