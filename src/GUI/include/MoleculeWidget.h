//
// Created by Michael Heuer on 12.11.17.
//

#ifndef AMOLQCPP_MOLECULEWIDGET_H
#define AMOLQCPP_MOLECULEWIDGET_H

#include <QWidget>
#include <Qt3DCore>
#include <Qt3DExtras>
#include <QVBoxLayout>
#include <QLabel>

class MoleculeWidget : public QWidget{
    Q_OBJECT
public:
    explicit MoleculeWidget(QWidget *parent = nullptr);
    Qt3DCore::QEntity* getRoot();

private:
    QVBoxLayout *layout_;
    Qt3DExtras::Qt3DWindow *qt3DWindow_;
    Qt3DCore::QEntity *root_;
    Qt3DExtras::QOrbitCameraController *cameraController_;
    QLabel* infoText_;
};

#endif //AMOLQCPP_MOLECULEWIDGET_H
