/* Copyright (C) 2018 Leonard Reuter.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include "Serialization.h"

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


    ElectronsVector electrons = YAML::LoadFile("LR_H2_artificial_start.json")["Electrons"].as<ElectronsVector>();

    double energy = -1.1745;

    LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>
            lagrangeProblem(electronicWaveFunctionProblem,potentialProblem,energy);

    GradientSqMagnitudeProblem<LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>> mainprob(lagrangeProblem);

    double lambdaInit = 0;

    Eigen::VectorXd x(electrons.positionsVector().asEigenVector());

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