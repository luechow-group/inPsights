//
// Created by Leonard Reuter on 14.03.18.
//

#ifndef AMOLQCPP_VISUALIZATION_H
#define AMOLQCPP_VISUALIZATION_H


#include <Qt3DCore>
#include <Qt3DRender>
#include <Qt3DExtras>

#include <QtWidgets/QApplication>

#include "MoleculeWidget.h"
#include "AtomsVector3D.h"
#include "ElectronsVector3D.h"

#include "ParticlesVectorPath3D.h"


namespace Visualization{
    int visualizeOptPath(int &argc, char **argv,
                         const AtomsVector &atoms,
                         const ElectronsVectorCollection &optimizationPath,
                         const unsigned long &nwanted = 300) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC,"C");


        ElectronsVectorCollection shortenedPath(optimizationPath[0]);

        double optPathLength = 0.0;
        Eigen::VectorXd pathLengthVector =
                Eigen::VectorXd::Zero(optimizationPath.numberOfEntities());

        for (unsigned long i = 1; i < optimizationPath.numberOfEntities(); i++){
            optPathLength += optimizationPath.norm(i,i - 1);
            pathLengthVector[i] = optPathLength;
        }

        double stepLength = optPathLength / nwanted;
        ElectronsVector elecVector;
        unsigned long index = 0;
        for (unsigned long i = 0; i <= nwanted; i++) {
            index = 0;
            for (unsigned long j = 1; j < optimizationPath.numberOfEntities(); j++){
                if (fabs(pathLengthVector[j] - i * stepLength) <
                    fabs(pathLengthVector[index] - i * stepLength)){
                    index = j;
                }
            }
            shortenedPath.append(optimizationPath[index]);
        }

        // Visualization
        MoleculeWidget moleculeWidget;
        Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

        AtomsVector3D(root, atoms);

        // Plot the starting point
        ElectronsVector3D(root, optimizationPath[-1], false);

        // Plot the optimization path
        ParticlesVectorPath3D(root, shortenedPath);

        return app.exec();
    }
}

#endif //AMOLQCPP_VISUALIZATION_H
