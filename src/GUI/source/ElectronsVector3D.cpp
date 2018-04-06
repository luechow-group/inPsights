//
// Created by Michael Heuer on 12.11.17.
//

#include "ElectronsVector3D.h"
#include "Electron3D.h"
#include "Bond3D.h"

#include <QPhongMaterial>
#include <QExtrudedTextMesh>
#include <AtomsVector.h>

ElectronsVector3D::ElectronsVector3D(Qt3DCore::QEntity *root,
                                     const ElectronsVector &electronsVector,
                                     bool showIndicesQ) {
    drawElectrons(root, electronsVector, showIndicesQ);
}

ElectronsVector3D::ElectronsVector3D(Qt3DCore::QEntity *root, const AtomsVector &atomsVector,
                                     const ElectronsVector &electonsVector, bool showIndicesQ)
        : ElectronsVector3D(root, electonsVector, showIndicesQ)
{
    drawConnections(root, atomsVector, electonsVector);
}

void ElectronsVector3D::drawElectrons(Qt3DCore::QEntity *root,
                                      const ElectronsVector &electronsVector,
                                      bool showIndicesQ) {

    // Draw electrons
    std::vector<Electron3D> electrons3D;

    for (long i = 0; i < electronsVector.numberOfEntities(); ++i) {
        Eigen::Vector3d vec = electronsVector[i].position();
        auto qvector3d = QVector3D(float(vec[0]), float(vec[1]), float(vec[2]));
        electrons3D.emplace_back(Electron3D(root, qvector3d,
                                            electronsVector.spinTypesVector()[i]));

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

void ElectronsVector3D::drawConnections(Qt3DCore::QEntity *root,
                                        const AtomsVector &atomsVector,
                                        const ElectronsVector &electronsVector) {


    std::vector<long> electronsNotInNuclei;

    for (long i = 0; i < electronsVector.numberOfEntities(); ++i) {
        bool notAtNuclei = true;
        for (int k = 0; k < atomsVector.numberOfEntities(); ++k) {
            if (Particle::distance(electronsVector[i], atomsVector[k]) < 0.1) {
                notAtNuclei *= false;
            }
        }
        if (notAtNuclei) electronsNotInNuclei.push_back(i);
    }


    for (long i = 0; i < electronsNotInNuclei.size(); ++i) {
        for (long j = i+1; j < electronsNotInNuclei.size(); ++j) {

            auto e1 = electronsVector[electronsNotInNuclei[i]];
            auto e2 = electronsVector[electronsNotInNuclei[j]];

            if (Particle::distance(e1, e2) < 2) {
                auto p1 = e1.position();
                auto p2 = e2.position();
                auto q1 = QVector3D(p1[0], p1[1], p1[2]);
                auto q2 = QVector3D(p2[0], p2[1], p2[2]);

                if (e1.spinType() == e2.spinType())
                    Cylinder(root, Qt::magenta, {q1,q2}, 0.005,0.5);
                else if (e1.spinType() != e2.spinType())
                    Cylinder(root, Qt::green, {q1,q2}, 0.005,0.5);
                else
                    Cylinder(root, Qt::cyan, {q1,q2}, 0.005, 0.25);
            }
        }
    }
}

