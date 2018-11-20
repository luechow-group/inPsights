//
// Created by Michael Heuer on 12.11.17.
//

#ifndef INPSIGHTS_MOLECULEWIDGET_H
#define INPSIGHTS_MOLECULEWIDGET_H

#include <QWidget>
#include <Qt3DCore>
#include <Qt3DExtras>
#include <QVBoxLayout>
#include <QLabel>

#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>
#include <Line3D.h>
class MoleculeWidget : public QWidget{
    Q_OBJECT
public:
    explicit MoleculeWidget(QWidget *parent = nullptr);
    Qt3DCore::QEntity* getRoot();

    void setMolecule(
            const AtomsVector& atoms,
            const std::pair<std::vector<ElectronsVector>,YAML::Node>& clusterData,
            bool drawConnections,
            bool drawSpinCorrelations,
            double spinCorrelationThreshold){

        moleculeEntity_->deleteLater();
        moleculeEntity_ = new Qt3DCore::QEntity(root_);

        atomsVector3D_ = new AtomsVector3D(moleculeEntity_, atoms);

        auto electrons = clusterData.first[0]; // Plot all?

        if(drawConnections)
            electronsVector3D_ = new ElectronsVector3D(moleculeEntity_, atoms, electrons);
        else
            electronsVector3D_ = new ElectronsVector3D(moleculeEntity_, electrons);

        if(drawSpinCorrelations) {
            for (int i = 0; i < electrons.numberOfEntities(); ++i) {
                for (int j = i + 1; j < electrons.numberOfEntities(); ++j) {

                    auto corr = clusterData.second[i][j][0].as<double>();
                    if (std::abs(corr) >= spinCorrelationThreshold) {

                        QColor color;
                        if(corr > 0)
                            color = QColor::fromRgb(255,0,255);
                        else
                            color = QColor::fromRgb(0,255,0);

                        QVector3D start, end;
                        start.setX(electrons.positionsVector()[i].x()); //TODO use helper
                        start.setY(electrons.positionsVector()[i].y());
                        start.setZ(electrons.positionsVector()[i].z());
                        end.setX(electrons.positionsVector()[j].x());
                        end.setY(electrons.positionsVector()[j].y());
                        end.setZ(electrons.positionsVector()[j].z());

                        new Line3D(moleculeEntity_, color, {start, end}, std::abs(corr));
                    }
                }
            }
        }
    }

private:
    QVBoxLayout *layout_;
    Qt3DExtras::Qt3DWindow *qt3DWindow_;
    Qt3DCore::QEntity *root_, *moleculeEntity_;
    Qt3DExtras::QOrbitCameraController *cameraController_;
    QLabel* infoText_;


    AtomsVector3D* atomsVector3D_;
    ElectronsVector3D* electronsVector3D_;
};

#endif //INPSIGHTS_MOLECULEWIDGET_H
