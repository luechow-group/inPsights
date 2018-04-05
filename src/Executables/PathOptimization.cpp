//
// Created by Michael Heuer on 07.11.17.
//

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

#include "Visualization.h"

bool handleCommandlineArguments(int argc, char **argv,
                                std::string &wavefunctionFilename,
                                std::string &electronsVectorFilename,
                                bool &showGui) {
    if (argc < 3) {
        std::cout << "Usage: \n"
                  << "Argument 1: wavefunction filename (.wf)\n"
                  << "Argument 2: electron collection filename (.json)\n"
                  << "Argument 3 (Optional): display the gui (gui)" << std::endl;
        std::cout << "Ethylene-em-5.wf LD_Ethlyen_Start.json gui" << std::endl;
        return false;
    } else if (argc >= 3) {
        wavefunctionFilename = argv[1];
        electronsVectorFilename = argv[2];
        if (argc > 3) showGui = (std::string(argv[3]) == "gui");
        return true;
    }
}

int main(int argc, char *argv[]) {
    std::string wavefunctionFilename; //= "BH3_Exp-em.wf"; // overwrite command line
    std::string electronsVectorFilename; //= "BH3_Max0.json"; // overwrite command line
    bool showGui = true;

    if( wavefunctionFilename.empty() && electronsVectorFilename.empty()) {
        bool inputArgumentsFoundQ =
                handleCommandlineArguments(argc, argv, wavefunctionFilename, electronsVectorFilename, showGui);
        if(!inputArgumentsFoundQ) return 0;
    }

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wavefunctionFilename);
    CollectionParser collectionParser;
    auto ac = electronicWaveFunctionProblem.getAtomsVector();
    auto ec = collectionParser.electronsVectorFromJson(electronsVectorFilename);
    std::cout << ac << std::endl;
    std::cout << ec << std::endl;

    Eigen::VectorXd x(ec.positionsVector().positionsAsEigenVector());
    std::cout << x.transpose() << std::endl;
    Eigen::VectorXd grad(ec.numberOfEntities());
    electronicWaveFunctionProblem.putElectronsIntoNuclei(x,grad);


    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();

    cppoptlib::BfgsUmrigarSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    crit.gradNorm = 1e-6;
    crit.iterations = 1000;
    solver.setStopCriteria(crit);
    //solver.setMaxStepLength(1.0);
    //solver.setSteepestDescentRate(0.1);
    //solver.setDistanceCriteriaUmrigar(0.1);

    solver.minimize(electronicWaveFunctionProblem, x);
    std::cout << ElectronsVector(x,ec.spinTypesVector().spinTypesAsEigenVector())<<std::endl;


    if(showGui) {
        auto optimizationPath = electronicWaveFunctionProblem.getOptimizationPath();

        //prepend starting point
        optimizationPath.prepend(ec);

        //append end point
        optimizationPath.append(ElectronsVector(PositionsVector(x), optimizationPath.spinTypesVector()));

        return Visualization::visualizeOptPath(argc, argv, ac, optimizationPath);
    }
    return 0;
};
