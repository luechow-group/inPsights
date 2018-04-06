//
// Created by Michael Heuer on 07.11.17.
//

#include <CollectionParser.h>
#include <AmolqcFileImport/RefFileImporter.h>
#include <ElectronicWaveFunctionProblem.h>
#include <MoleculeWidget.h>
#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>

#include <QApplication>
#include <iostream>
#include <exception>

bool handleCommandlineArguments(int argc, char **argv,
                                std::string &filename, unsigned &k, unsigned &m) {
    if (argc < 4) {
        std::cout << "Usage: \n"
                  << "Argument 1: reference filename (.ref)" << std::endl
                  << "Argument 2: k number" << std::endl
                  << "Argument 3: m number" << std::endl
                  << "max.ref 1 2" << std::endl;
        return false;
    } else if (argc == 4) {
        filename = argv[1];
        int kint = std::atoi(argv[2]);
        int mint = std::atoi(argv[3]);
        try {
            if (".ref" != filename.substr( filename.length() - 4 ))
                throw std::invalid_argument("Invalid file type '" + filename + "'.");
            if (kint < 1)
                throw std::invalid_argument("The k index must be a positive integer but is " + std::string(argv[2]));
            if (mint < 1)
                throw std::invalid_argument("The m index must be a positive integer but is " + std::string(argv[3]));
        }
        catch (std::invalid_argument &e){
            std::cout << e.what() << std::endl;
            abort();
        }

        k = unsigned(kint);
        m = unsigned(mint);
        return true;
    }
}

int main(int argc, char *argv[]) {
    std::string refFilename;
    unsigned k,m;

    if(!handleCommandlineArguments(argc, argv, refFilename,k,m)) return 0;

    RefFileImporter refFileImporter(refFilename);
    auto atoms = refFileImporter.getAtomsVector();
    auto electrons = refFileImporter.getMaximaStructure(k,m);

    // Visualization
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    AtomsVector3D(root, atoms);
    ElectronsVector3D(root, atoms, electrons, false);

    return app.exec();
};
