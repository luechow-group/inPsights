//
// Created by heuer on 10.05.17.
//

#ifndef AMOLQCGUI_ELECTRON_H
#define AMOLQCGUI_ELECTRON_H

#include "Sphere.h"
#include "ElementInfo.h"

namespace Spin {
    enum SpinType { Alpha, Beta };

    static QColor QColorFromSpinType(const SpinType& spinType){
      switch (spinType){
        case Alpha:
          return Qt::red;
        case Beta:
          return Qt::blue;
      }
    }
}

class Electron : public Sphere {
    //Q_OBJECT
public:
    Electron(const Electron& electron);
    Electron(Qt3DCore::QEntity *root,
         const QVector3D& location,
         const Spin::SpinType& spinType);

    Spin::SpinType getSpinType() const { return spinType_; };

    Qt3DRender::QObjectPicker *picker;

//public slots:
    //void onPressed(bool pressed);

private:
    const Spin::SpinType spinType_;
};

#endif //AMOLQCGUI_ELECTRON_H
