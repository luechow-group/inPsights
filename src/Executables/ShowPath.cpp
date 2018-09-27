//
// Created by Michael Heuer on 27.09.18.
//

#include <iostream>
#include <AmolqcFileImport/WfFileImporter.h>
#include <AmolqcFileImport/OptimizationPathFileImporter.h>
#include <Visualization.h>

bool handleCommandlineArguments(int argc, char **argv,
                                std::string &pathFilename,
                                unsigned & pathId) {
    if (argc < 2) {
        std::cout << "Usage: \n"
                  << "Argument 1: path filename (.300)\n"
                  << "Argument 2: path id"
                  << std::endl;
        std::cout << "Ethane.300 1" << std::endl;
        return false;
    } else if (argc == 3) {
        pathFilename = argv[1];
        int kint = std::atoi(argv[2]);
        try {
            if (".300" != pathFilename .substr( pathFilename .length() - 4 ))
                throw std::invalid_argument("Invalid path file type '" + pathFilename  + "'.");
            if (kint < 1)
                throw std::invalid_argument("The path index must be a positive integer but is " + std::string(argv[2]));
        }
        catch (std::invalid_argument &e){
            std::cout << e.what() << std::endl;
            abort();
        }

        pathId = unsigned(kint);
        return true;
    } else {
        throw std::invalid_argument("Too many arguments");
    }
};

int main(int argc, char *argv[]) {
    std::string pathFilename;
    unsigned pathId;

    if (pathFilename.empty()) {
        bool inputArgumentsFoundQ =
                handleCommandlineArguments(argc, argv, pathFilename, pathId);
        if (!inputArgumentsFoundQ) return 0;
    }


    OptimizationPathFileImporter optimizationPathFileImporter(pathFilename,1);//TODO REMOVE MULTIPLICITY
    auto path = optimizationPathFileImporter.getPath(pathId);
    auto atoms = optimizationPathFileImporter.getAtomsVector();

    return Visualization::visualizeOptPath(argc, argv, atoms, path);
}
