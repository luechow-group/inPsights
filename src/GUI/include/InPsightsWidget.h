#include <utility>

//
// Created by Michael Heuer on 18.11.18.
//

#ifndef INPSIGHTS_INPSIGHTSWIDGET_H
#define INPSIGHTS_INPSIGHTSWIDGET_H

#include <QWidget>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QSpinBox>
#include <MoleculeWidget.h>
#include <QTimer>
#include <QSplashScreen>
#include <QCheckBox>

class InPsightsWidget : public QWidget{
    Q_OBJECT
public:

    explicit InPsightsWidget(
            AtomsVector atomsVector,
            std::vector<ElectronsVector> electronsVectorCollection,
            QWidget *parent = nullptr
            )
    :
    QWidget(parent),
    moleculeWidget_(new MoleculeWidget(this)),
    idSpinBox_(new QSpinBox(this)),
    spinConnectionsCheckBox_(new QCheckBox("Spin Connections",this)),
    atomsVector_(std::move(atomsVector)),
    electronsVectorCollection_(std::move(electronsVectorCollection))
    {
        auto splashScreen = createSplashScreen();

        auto hbox = new QHBoxLayout(this);
        auto vbox = new QVBoxLayout(this);
        auto gbox = new QGroupBox("Settings:",this);

        setLayout(hbox);
        resize(800,800);
        hbox->addWidget(moleculeWidget_, Qt::AlignLeft);
        hbox->addWidget(gbox);
        gbox->setLayout(vbox);
        vbox->addWidget(idSpinBox_);
        vbox->addWidget(spinConnectionsCheckBox_);


        moleculeWidget_->setMolecule(atomsVector_,electronsVectorCollection_[0]);

        QObject::connect(idSpinBox_, SIGNAL(valueChanged(int)),
                         this, SLOT(redrawMolecule(int)));
        QObject::connect(spinConnectionsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(onConnectionsStateChanged(int)));

        idSpinBox_->setRange(0, int(electronsVectorCollection_.size()-1));
        idSpinBox_->setSingleStep(1);
        idSpinBox_->setValue(0);

        QTimer::singleShot(2000, splashScreen, SLOT(close()));
        QTimer::singleShot(2000, this, SLOT(show()));
    }

public slots:
    void redrawMolecule(int id){
        moleculeWidget_->setMolecule(atomsVector_, electronsVectorCollection_[id],
                                     spinConnectionsCheckBox_->checkState() == Qt::CheckState::Checked);
    };

    void onConnectionsStateChanged(int state){
        moleculeWidget_->setMolecule(atomsVector_, electronsVectorCollection_[idSpinBox_->value()],
                                     spinConnectionsCheckBox_->checkState() == Qt::CheckState::Checked);
    };


private:
    MoleculeWidget *moleculeWidget_;
    QSpinBox *idSpinBox_;
    QCheckBox *spinConnectionsCheckBox_;

    AtomsVector atomsVector_;
    std::vector<ElectronsVector> electronsVectorCollection_;



    QSplashScreen* createSplashScreen(){
        auto splashScreen = new QSplashScreen();

        QImage file(":inPsights.png");
        QPixmap qPixmap;
        qPixmap.convertFromImage(file);
        splashScreen->setPixmap(qPixmap);
        splashScreen->show();
        splashScreen->showMessage("Version 1.0.0", Qt::AlignRight, Qt::lightGray);

        return splashScreen;
    }
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
