//
// Created by heuer on 10.05.17.
//

#ifndef INPSIGHTS_ELECTRON3D_H
#define INPSIGHTS_ELECTRON3D_H

#include "Sphere.h"
#include "ElementInfo.h"
#include "SpinType.h"

class Electron3D : public Sphere {
    //Q_OBJECT
public:
    Electron3D(const Electron3D& electron);
    Electron3D(Qt3DCore::QEntity *root,
         const QVector3D& location,
         const Spin& spinType);

    Spin getSpinType() const { return spinType_; };

    Qt3DRender::QObjectPicker *picker;

//public slots:
    //void onPressed(bool pressed);

private:
    const Spin spinType_;
};

#endif //INPSIGHTS_ELECTRON3D_H
