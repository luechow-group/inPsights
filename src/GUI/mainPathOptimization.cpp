//
// Created by Michael Heuer on 07.11.17.
//

#include <QApplication>
#include <iostream>

#include "OptimizationPathFileImporter.h"
#include "WfFileImporter.h"

#include "AtomCollection3D.h"
#include "ElectronCollection3D.h"
#include "MoleculeWidget.h"


#include "ParticleCollectionPath3D.h"

#include <Qt3DRender>
#include <Qt3DExtras>
#include <QtWidgets/QApplication>
#include <QtWidgets>




int main(int argc, char *argv[]) {

    QApplication app(argc, argv);


    std::string filename = "Diborane.wf";
    WfFileImporter wfFileImporter(filename);
    auto ac = wfFileImporter.getAtomCollection();

    OptimizationPathFileImporter optimizationPathFileImporter("Diborane-Paths.300",1);
    auto ecs = optimizationPathFileImporter.getPath(1);

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    AtomCollection3D(root, ac);
    ElectronCollection3D(root, ecs.getElectronCollection(0));

    ParticleCollectionPath3D(root, ecs);

    return app.exec();
};