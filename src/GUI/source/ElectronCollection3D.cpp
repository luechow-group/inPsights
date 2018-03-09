//
// Created by Michael Heuer on 12.11.17.
//

#include "ElectronCollection3D.h"
#include "Electron3D.h"
#include "Bond3D.h"

#include <QPhongMaterial>
#include <QExtrudedTextMesh>

ElectronCollection3D::ElectronCollection3D(Qt3DCore::QEntity *root, const ElectronCollection &electronCollection,
                                           bool showIndicesQ) {

    std::vector<Electron3D> electrons3D;

    // Draw electrons
    for (long i = 0; i < electronCollection.numberOfEntities(); ++i) {
        Eigen::Vector3d vec= electronCollection[i].position();
        auto qvector3d = QVector3D(float(vec[0]),float(vec[1]),float(vec[2]));
        electrons3D.emplace_back(Electron3D(root, qvector3d,
                                            electronCollection.spinTypeCollection()[i]));

        // Draw Text
        if(showIndicesQ){
            auto *textMaterial = new Qt3DExtras::QPhongMaterial(root);

            auto *text = new Qt3DCore::QEntity(root);
            auto *textMesh = new Qt3DExtras::QExtrudedTextMesh();

            auto *textTransform = new Qt3DCore::QTransform();
            QFont font(QString("Arial"), 12, 0, false);
            QVector3D shift;

            if(electrons3D[i].getSpinType() == Spin::SpinType::alpha) shift = QVector3D(0.0f,0.07f,0.07f);
            else shift = QVector3D(-0.0f,-0.07f,-0.07f);

            textTransform->setTranslation(qvector3d+shift);
            textTransform->setRotationY(90);
            textTransform->setScale(0.1f);
            textMesh->setDepth(0.1f);
            textMesh->setFont(font);
            textMesh->setText(QString::fromStdString(std::to_string(i+1)));
            textMaterial->setDiffuse(electrons3D[i].getColor());

            text->addComponent(textMaterial);
            text->addComponent(textMesh);
            text->addComponent(textTransform);
        }
    }

}




