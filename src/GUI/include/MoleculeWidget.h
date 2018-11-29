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

#include <Metrics.h>
#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>
#include <Line3D.h>
#include <Bond3D.h>
#include <Cylinder.h>
#include <Electron3D.h>
#include <ClusterData.h>
#include <AtomsVectorLinkedElectronsVector.h>


class MoleculeWidget : public QWidget{
    Q_OBJECT
public:
    explicit MoleculeWidget(QWidget *parent = nullptr);

    Qt3DCore::QEntity* getRoot();

    //TODO implement
    void deleteConnections() {};

    void drawConnections() { //TODO How to make them deletable?
        for (auto &mapItem : activeElectronsVectorsMap_) {
            auto electronsVector3D = mapItem.second;
            AtomsVectorLinkedElectronsVector linkedElectronsVector(
                    sharedAtomsVector_,
                    electronsVector3D->electronsVector_);

            double coreThreshold = 0.1;
            double maxDistance = 1.6;

            auto valenceElectronIndices = linkedElectronsVector.valenceElectronsIndices(coreThreshold);

            for (auto it = valenceElectronIndices.begin(); it != valenceElectronIndices.end(); ++it) {
                for (auto jt = it + 1; jt != valenceElectronIndices.end(); ++jt) {
                    auto e1 = linkedElectronsVector[*it];
                    auto e2 = linkedElectronsVector[*jt];

                    auto q1 = GuiHelper::toQVector3D(e1.position());
                    auto q2 = GuiHelper::toQVector3D(e2.position());

                    if (Metrics::distance(e1.position(), e2.position()) < maxDistance) {
                        if (e1.type() == e2.type())
                            if (e1.type() == Spin::alpha)
                                Cylinder(electronsVector3D, GuiHelper::QColorFromType<Spin>(Spin::alpha), {q1, q2}, 0.015,
                                         0.5);
                            else if (e1.type() == Spin::beta)
                                Cylinder(electronsVector3D, GuiHelper::QColorFromType<Spin>(Spin::beta), {q1, q2}, 0.015,
                                         0.5);
                            else
                                Cylinder(electronsVector3D, GuiHelper::QColorFromType<Spin>(Spin::none), {q1, q2}, 0.015,
                                         0.5);
                        else
                            Line3D(electronsVector3D, Qt::black, {q1, q2}, 0.25);
                    }
                }
            }
        }
    }

    // TODO should be individually selectable for each list item
    /*void drawSpinCorrelations() {
        for ( auto it = activeElectronsVectorsMap_.begin(); it != activeElectronsVectorsMap_.end(); it++ ){
            if(drawSpinCorrelations) {
                for (int i = 0; i < electrons.numberOfEntities(); ++i) {
                    for (int j = i + 1; j < electrons.numberOfEntities(); ++j) {

                        auto corr = clusterData.SeeStats_.mean()(i,j);
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
    }*/

    void drawAtoms(){
        atomsVector3D_ = new AtomsVector3D(moleculeEntity_, *sharedAtomsVector_);
    }

    void drawBonds() {
        // Draw bonds
        auto bondDrawingLimit = float(1.40 * 1e-10 / AU::length);

        for (long i = 0; i < sharedAtomsVector_->numberOfEntities(); ++i) {
            for (long j = i + 1; j < sharedAtomsVector_->numberOfEntities(); ++j) {
                auto a1 = sharedAtomsVector_->operator[](i);
                auto a2 = sharedAtomsVector_->operator[](j);
                auto distance = Metrics::distance(a1.position(), a2.position())
                        - (Elements::ElementInfo::vdwRadius(a1.type()) + Elements::ElementInfo::vdwRadius(a2.type()))/10.0;

                Elements::ElementInfo::vdwRadius(a1.type());
                if ( distance < bondDrawingLimit) {
                    new Bond3D(moleculeEntity_, a1, a2);
                }
            }
        }
    }

    void addElectronsVector(const ElectronsVector& electronsVector, int id = 0){
        activeElectronsVectorsMap_.emplace(id,new ElectronsVector3D(moleculeEntity_, electronsVector));
    }

    void removeElectronsVector(int id = 0){
        activeElectronsVectorsMap_[id]->deleteLater();
        activeElectronsVectorsMap_.erase(id);
    }

    void setSharedAtomsVector(AtomsVector atomsVector) {
        sharedAtomsVector_ = std::make_shared<AtomsVector >(std::move(atomsVector));
    }

private:
    QVBoxLayout *layout_;
    Qt3DExtras::Qt3DWindow *qt3DWindow_;
    Qt3DCore::QEntity *root_, *moleculeEntity_;
    Qt3DExtras::QOrbitCameraController *cameraController_;
    QLabel* infoText_;

    std::shared_ptr<AtomsVector> sharedAtomsVector_;
    AtomsVector3D* atomsVector3D_;
    std::map<int, ElectronsVector3D*> activeElectronsVectorsMap_;
};

#endif //INPSIGHTS_MOLECULEWIDGET_H
