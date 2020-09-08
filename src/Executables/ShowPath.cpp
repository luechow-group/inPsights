// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <iostream>
#include <WfFileImporter.h>
#include <OptimizationPathFileImporter.h>
#include <Visualization.h>
#include <MoleculeWidget.h>
#include <ParticlesVectorPath3D.h>

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

    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    OptimizationPathFileImporter optimizationPathFileImporter(pathFilename);
    auto path = optimizationPathFileImporter.getPath(pathId);
    auto atoms = optimizationPathFileImporter.getAtomsVector();

    long nWanted = 1000;
    ElectronsVectorCollection visualizationPath;
    if (nWanted < path.numberOfEntities())
        visualizationPath = Visualization::shortenPath(path, nWanted);
    else 
        visualizationPath = path;

    auto moleculeWidget = new MoleculeWidget();
    moleculeWidget->setSharedAtomsVector(atoms);
    moleculeWidget->addElectronsVector(visualizationPath[-1]); // Last ElectronsVector of the path
    moleculeWidget->drawAtoms();
    moleculeWidget->drawBonds();

    ParticlesVectorPath3D(moleculeWidget->getMoleculeEntity(), visualizationPath);

    moleculeWidget->fileInfoText_->setText(QString::fromStdString(pathFilename));
    moleculeWidget->resize(1024,768);
    moleculeWidget->show();

    return app.exec();
}
