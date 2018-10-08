//
// Created by heuer on 10.05.17.
//

#ifndef AMOLQCPP_ELECTRON3D_H
#define AMOLQCPP_ELECTRON3D_H

#include "Sphere.h"
#include "ElementInfo.h"
#include "SpinType.h"

// add QColorFromSpinType method to the namespace Spin declared in BaseLib
namespace Spins {
    static QColor QColorFromSpinType(const SpinType& spinType){
        switch (spinType){
            case SpinType::alpha:
                return Qt::red;
            case SpinType::beta:
                return Qt::blue;
            case SpinType::none:
                return Qt::darkMagenta;
            default:
                return Qt::darkMagenta;
        }
    }
}

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

#endif //AMOLQCPP_ELECTRON3D_H
