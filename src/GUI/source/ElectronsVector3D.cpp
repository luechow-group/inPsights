//
// Created by Michael Heuer on 12.11.17.
//

#include <ElectronsVector3D.h>
#include <Electron3D.h>
#include <Bond3D.h>
#include <Line3D.h>
#include <Metrics.h>
#include <QPhongMaterial>
#include <QExtrudedTextMesh>
#include <ParticlesVector.h>

ElectronsVector3D::ElectronsVector3D(Qt3DCore::QEntity *root,
                                     const ElectronsVector &electronsVector,
                                     bool showIndicesQ)
                                     : QEntity(root), electronsVector_(electronsVector) {
    drawElectrons(showIndicesQ);
}

ElectronsVector3D::ElectronsVector3D(Qt3DCore::QEntity *root, const AtomsVector &atomsVector,
                                     const ElectronsVector &electonsVector, bool showIndicesQ)
        : ElectronsVector3D(root, electonsVector, showIndicesQ)
{
    drawConnections(atomsVector, electonsVector);
}

void ElectronsVector3D::drawElectrons(bool showIndicesQ) {
    // Draw electrons
    std::vector<Electron3D*> electrons3D;

    for (long i = 0; i < electronsVector_.numberOfEntities(); ++i) {
        auto qvector3d = GuiHelper::toQVector3D(electronsVector_[i].position());
        electrons3D.emplace_back(new Electron3D(this, qvector3d, electronsVector_.typesVector()[i]));

        // Draw Text
        if(showIndicesQ){
            auto *textMaterial = new Qt3DExtras::QPhongMaterial(this);

            auto *text = new Qt3DCore::QEntity(this);
            auto *textMesh = new Qt3DExtras::QExtrudedTextMesh();

            auto *textTransform = new Qt3DCore::QTransform();
            QFont font(QString("Arial"), 12, 0, false);
            QVector3D shift;

            if(electrons3D[i]->getSpinType() == Spin::alpha) shift = QVector3D(0.0f,0.15f,0.15f);
            else shift = QVector3D(-0.0f,-0.15f,-0.15f);

            textTransform->setTranslation(qvector3d+shift);
            textTransform->setRotationY(90);
            textTransform->setScale(0.2f);
            textMesh->setDepth(0.1f);
            textMesh->setFont(font);
            textMesh->setText(QString::fromStdString(std::to_string(i+1)));
            textMaterial->setDiffuse(electrons3D[i]->getColor());

            text->addComponent(textMaterial);
            text->addComponent(textMesh);
            text->addComponent(textTransform);
        }
    }
}
//TODO DELETE
void ElectronsVector3D::drawConnections(const AtomsVector &atomsVector,
                                        const ElectronsVector &electronsVector) {
    std::vector<long> electronsNotInNuclei;

    for (long i = 0; i < electronsVector.numberOfEntities(); ++i) {
        bool notAtNuclei = true;
        for (long k = 0; k < atomsVector.numberOfEntities(); ++k) {
            if (Metrics::distance(electronsVector[i].position(), atomsVector[k].position()) < 0.1) {
                notAtNuclei = false;
            }
        }
        if (notAtNuclei) electronsNotInNuclei.push_back(i);
    }

    for (size_t i = 0; i < electronsNotInNuclei.size(); ++i) {
        for (size_t j = i+1; j < electronsNotInNuclei.size(); ++j) {

            auto e1 = electronsVector[electronsNotInNuclei[i]];
            auto e2 = electronsVector[electronsNotInNuclei[j]];

            double distanceThreshold = 1.6;
            if (Metrics::distance(e1.position(), e2.position()) < distanceThreshold) {
                auto p1 = e1.position();
                auto p2 = e2.position();
                auto q1 = QVector3D(p1[0], p1[1], p1[2]);
                auto q2 = QVector3D(p2[0], p2[1], p2[2]);

                if (e1.type() == e2.type())
                    if(e1.type() == Spin::alpha)
                        Cylinder(this,GuiHelper::QColorFromType<Spin>(Spin::alpha), {q1,q2}, 0.015,0.5);
                    else if(e1.type() == Spin::beta)
                        Cylinder(this,GuiHelper::QColorFromType<Spin>(Spin::beta), {q1,q2}, 0.015,0.5);
                    else
                        Cylinder(this,GuiHelper::QColorFromType<Spin>(Spin::none), {q1,q2}, 0.015,0.5);
                else
                    Line3D(this, Qt::black, {q1,q2}, 0.25);
            }
        }
    }
}

