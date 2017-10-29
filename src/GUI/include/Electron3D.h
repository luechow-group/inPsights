//
// Created by heuer on 10.05.17.
//

#ifndef AMOLQCGUI_ELECTRON3D_H
#define AMOLQCGUI_ELECTRON3D_H

#include "Sphere.h"
#include "ElementInfo.h"
#include "SpinType.h"

// add QColorFromSpinType method to the namespace Spin declared in BaseLib
namespace Spin {
    static QColor QColorFromSpinType(const SpinType& spinType){
      switch (spinType){
        case Alpha:
          return Qt::red;
        case Beta:
          return Qt::blue;
        case None:
          return Qt::black;
      }
    }
}

class Electron3D : public Sphere {
    //Q_OBJECT
public:
    Electron3D(const Electron3D& electron);
    Electron3D(Qt3DCore::QEntity *root,
         const QVector3D& location,
         const Spin::SpinType& spinType);

    Spin::SpinType getSpinType() const { return spinType_; };

    Qt3DRender::QObjectPicker *picker;

//public slots:
    //void onPressed(bool pressed);

private:
    const Spin::SpinType spinType_;
};

#endif //AMOLQCGUI_ELECTRON3D_H
