//
// Created by Leonard Reuter on 15.03.18.
//

#include "Visualization.h"
namespace Visualization {
    
    namespace { // anonymous namespace for private functions
        ElectronsVectorCollection shortenPath(const ElectronsVectorCollection &optimizationPath, const unsigned long &nwanted) {
            ElectronsVectorCollection visualizationPath(optimizationPath[0]);

            double optPathLength = 0.0;
            Eigen::VectorXd pathLengthVector =
                    Eigen::VectorXd::Zero(optimizationPath.numberOfEntities());

            for (unsigned long i = 1; i < optimizationPath.numberOfEntities(); i++){
                optPathLength += optimizationPath.norm(i,i - 1);
                pathLengthVector[i] = optPathLength;
            }

            double stepLength = optPathLength / nwanted;

            unsigned long index = 0;
            for (unsigned long i = 0; i <= nwanted; i++) {
                index = 0;
                for (unsigned long j = 1; j < optimizationPath.numberOfEntities(); j++){
                    if (fabs(pathLengthVector[j] - i * stepLength) <
                            fabs(pathLengthVector[index] - i * stepLength)){
                        index = j;
                    }
                }
                visualizationPath.append(optimizationPath[index]);
            }
        }
    }

    int visualizeOptPath(int &argc, char **argv,
                       const AtomsVector &atoms,
                       const ElectronsVectorCollection &optimizationPath,
                       const unsigned long &nwanted) {

        QApplication app(argc, argv);
        setlocale(LC_NUMERIC,"C");

        ElectronsVectorCollection visualizationPath;
        if (nwanted < optimizationPath.numberOfEntities()){
            visualizationPath = shortenPath(optimizationPath, nwanted);
        }
        else {
            visualizationPath = optimizationPath;
        }

        // Visualization
        MoleculeWidget moleculeWidget;
        Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

        AtomsVector3D(root, atoms);

        // Plot the end point
        ElectronsVector3D(root, visualizationPath[-1], false);

        // Plot the optimization path
        ParticlesVectorPath3D(root, visualizationPath);

        return app.exec();
    }
}