// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "MoleculeWidget.h"
#include "ParticlesVector3D.h"

#include "ParticlesVectorPath3D.h"
#include "Polyline.h"
#include "Visualization.h"

namespace Visualization {

    ElectronsVectorCollection shortenPath(const ElectronsVectorCollection &optimizationPath, long nwanted) {
        ElectronsVectorCollection visualizationPath(optimizationPath[0]);

        double optPathLength = 0.0;
        Eigen::VectorXd pathLengthVector =
                Eigen::VectorXd::Zero(optimizationPath.numberOfEntities());

        for (long i = 1; i < optimizationPath.numberOfEntities(); i++){
            optPathLength += optimizationPath.norm(i,i - 1);
            pathLengthVector[i] = optPathLength;
        }

        double stepLength = optPathLength / (nwanted - 1);

        long index = 0;
        // start at 1 because visualization Path already contains optpath[0]
        for (long i = 1; i < nwanted; i++) {
            index = 0;
            for (long j = 1; j < optimizationPath.numberOfEntities(); j++){
                if (fabs(pathLengthVector[j] - i * stepLength) <
                        fabs(pathLengthVector[index] - i * stepLength)){
                    index = j;
                }
            }
            visualizationPath.append(optimizationPath[index]);
        }
        return visualizationPath;
    }

    int visualizeOptPath(int &argc, char **argv,
                       const AtomsVector &atoms,
                       const ElectronsVectorCollection &optimizationPath,
                       long nwanted) {

        QApplication app(argc, argv);
        setlocale(LC_NUMERIC,"C");

        ElectronsVectorCollection visualizationPath;
        if (long(nwanted) < optimizationPath.numberOfEntities()){
            visualizationPath = shortenPath(optimizationPath, nwanted);
        }
        else {
            visualizationPath = optimizationPath;
        }

        // Visualization
        auto moleculeWidget = new MoleculeWidget();
        Qt3DCore::QEntity *root = moleculeWidget->getMoleculeEntity();

        AtomsVector3D(root, atoms);

        // Plot the end point
        ElectronsVector3D(root, visualizationPath[-1]);

        // Plot the optimization path
        ParticlesVectorPath3D(root, visualizationPath);

        return app.exec();
    }
}

void Visualization::drawEigenVector(Qt3DCore::QEntity *root,
                                    const Eigen::MatrixXd &eigenvectors,
                                    const Eigen::VectorXd &origin,
                                    int eigenvectorIndex) {

    assert(eigenvectorIndex < eigenvectors.cols() && "Index must smaller than the number of eigenvectors.");
    assert(eigenvectorIndex >= 0 && "Index must be positive.");
    assert(eigenvectors.rows()%3 == 0 && "The number of rows must be dividable by 3.");
    assert(eigenvectors.cols()%3 == 0 && "The number of columns must be dividable by 3.");


    for (int i = 0; i < eigenvectors.rows()/3; ++i) {
        QVector3D v1(GuiHelper::toQVector3D(Eigen::Vector3d(origin.segment(i,3))));
        QVector3D v2 = v1;
        QVector3D ev(GuiHelper::toQVector3D(Eigen::Vector3d(eigenvectors.col(eigenvectorIndex).segment(i,3))));
        v2 += ev;
        std::vector<QVector3D> points = {v1, v2};
        Polyline pl(root, QColor(Qt::black), points, 0.01, true);
    }
}

