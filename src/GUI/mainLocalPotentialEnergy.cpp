#include "CollectionParser.h"
#include "WfFileImporter.h"

#include "solver/gradientdescentsolver.h"

#include "ElectronicWaveFunctionProblem.h"
#include "PotentialProblem.h"
#include "LagrangeProblem.h"
#include "GradientSqMagnitudeProblem.h"

#include "Visualization.h"

using namespace Eigen;

int main(int argc, char *argv[]) {

    bool showGui = true;

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem("H2.wf");
    AtomsVector nuclei = electronicWaveFunctionProblem.getAtomsVector();

    PotentialProblem potentialProblem(nuclei);

    CollectionParser collectionParser;

    //ElectronsVector electrons = collectionParser.electronsVectorFromJson("LD_Diboran_Start.json");

    ElectronsVector electrons;
    electrons.append(Electron(Vector3d(0,0,0.1),Spin::SpinType::alpha));
    electrons.append(Electron(Vector3d(0,0,-0.1),Spin::SpinType::beta));

    double energy = -1.7;

    LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>
            lagrangeProblem(electronicWaveFunctionProblem,potentialProblem,energy);

    GradientSqMagnitudeProblem<LagrangeProblem<ElectronicWaveFunctionProblem,PotentialProblem>> mainprob(lagrangeProblem);

    double lambdaInit = -1;

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
                                                optimizationPath.spinTypesVector()));

        return Visualization::visualizeOptPath(argc, argv, optimizationPath);

    }

    return 0;
}