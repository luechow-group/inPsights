//
// Created by Michael Heuer on 27.09.18.
//

#include <iostream>
#include <AmolqcFileImport/WfFileImporter.h>
#include <AmolqcFileImport/OptimizationPathFileImporter.h>
#include <Visualization.h>

bool handleCommandlineArguments(int argc, char **argv,
                                std::string &wavefunctionFilename,
                                std::string &pathFilename,
                                int & pathId) {
    if (argc < 3) {
        std::cout << "Usage: \n"
                  << "Argument 1: wavefunction filename (.wf)\n"
                  << "Argument 2: path filename (.300)\n"
                  << "Argument 3: path id"
                  << std::endl;
        std::cout << "Ethane.wf Ethane.300 " << std::endl;
        return false;
    } else if (argc >= 3) {
        wavefunctionFilename = argv[1];
        pathFilename = argv[2];
        pathId = std::atoi(argv[3]);
        return true;
    } else
        return false;
}

int main(int argc, char *argv[]) {
    std::string wavefunctionFilename;
    std::string pathFilename;
    int pathId;

    if (wavefunctionFilename.empty() && pathFilename.empty()) {
        bool inputArgumentsFoundQ =
                handleCommandlineArguments(argc, argv, wavefunctionFilename, pathFilename, pathId);
        if (!inputArgumentsFoundQ) return 0;
    }

    WfFileImporter wfFileImporter(wavefunctionFilename);
    auto atoms = wfFileImporter.getAtomsVector();

    auto multiplicity = wfFileImporter.getMultiplicity();

    OptimizationPathFileImporter optimizationPathFileImporter(pathFilename,multiplicity);
    auto path = optimizationPathFileImporter.getPath(pathId);

    return Visualization::visualizeOptPath(argc, argv, atoms, path);
}
