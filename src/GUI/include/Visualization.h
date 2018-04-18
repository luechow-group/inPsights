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
                         const Eigen::VectorXd& origin, int eigenvectorIndex);



}

#endif //AMOLQCPP_VISUALIZATION_H
