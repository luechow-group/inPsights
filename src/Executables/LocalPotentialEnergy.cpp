//
// Created by Leonard Reuter on 07.03.18.
//

#include "CollectionParser.h"

#include "solver/gradientdescentsolver.h"

#include "ElectronicWaveFunctionProblem.h"
#include "PotentialProblem.h"
#include "LagrangeProblem.h"
#include "GradientSqMagnitudeProblem.h"

#include "Visualization.h"

using namespace Eigen;

int main(int argc, char *argv[]) {

    bool showGui = true;

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("H2.wf",false);
    AtomsVector nuclei = electronicWaveFunctionProblem.getAtomsVector();

    PotentialProblem potentialProblem(nuclei);


    ElectronsVector electrons = CollectionParser::electronsVectorFromJson(
            CollectionParser::readJSON("LR_H2_artificial_start.json"));

    double energy = -1.1745;

    LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>
            lagrangeProblem(electronicWaveFunctionProblem,potentialProblem,energy);

    GradientSqMagnitudeProblem<LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>> mainprob(lagrangeProblem);

    double lambdaInit = 0;

    Eigen::VectorXd x(electrons.positionsVector().positionsAsEigenVector());

    Eigen::VectorXd y = x;
    y.conservativeResize(y.size()+1,Eigen::NoChange);
    y[y.size() - 1] = lambdaInit;

    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
 
    cppoptlib::GradientDescentSolver<GradientSqMagnitudeProblem<LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>>> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);

    solver.minimize(mainprob, y);

    if(showGui) {

        // Prepare the optimization path for visualization
        auto optimizationPath = mainprob.getProblem().getProblem().getOptimizationPath();

        //prepend starting position
        optimizationPath.prepend(electrons);

        //append end position
        optimizationPath.append(ElectronsVector(PositionsVector(y.head(y.size() - 1)),
                                                optimizationPath.typesVector()));

        return Visualization::visualizeOptPath(argc, argv, nuclei, optimizationPath);
    }

    return 0;
}