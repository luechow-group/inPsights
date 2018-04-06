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
#include "Polyline.h"


namespace Visualization{
    int visualizeOptPath(int &argc, char **argv,
                         const AtomsVector &atoms,
                         const ElectronsVectorCollection &optimizationPath,
                         const unsigned long &nwanted = 300);

    ElectronsVectorCollection shortenPath(const ElectronsVectorCollection &optimizationPath,
                                const unsigned long &nwanted);

    void drawEigenVector(Qt3DCore::QEntity *root,
                         const Eigen::MatrixXd eigenvectors,
                         const Eigen::VectorXd& origin, int eigenvectorIndex){

        assert(eigenvectorIndex < eigenvectors.cols() && "Index must smaller than the number of eigenvectors.");
        assert(eigenvectorIndex >= 0 && "Index must be positive.");
        assert(eigenvectors.rows()%3 == 0 && "The number of rows must be dividable by 3.");
        assert(eigenvectors.cols()%3 == 0 && "The number of columns must be dividable by 3.");


        for (int i = 0; i < eigenvectors.rows()/3; ++i) {
            QVector3D v1(origin(i*3+0),
                         origin(i*3+1),
                         origin(i*3+2));
            QVector3D v2 = v1;
            QVector3D ev(eigenvectors.col(eigenvectorIndex)(i*3+0),
                         eigenvectors.col(eigenvectorIndex)(i*3+1),
                         eigenvectors.col(eigenvectorIndex)(i*3+2));
            v2 += ev;
            std::vector<QVector3D> points = {v1, v2};
            Polyline pl(root, QColor(Qt::black), points, 0.01, true);
        }
    }



}

#endif //AMOLQCPP_VISUALIZATION_H
