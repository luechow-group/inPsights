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
#include <QSlider>
#include <QFileDialog>

class InPsightsWidget : public QWidget {
Q_OBJECT
public:

    explicit InPsightsWidget(QWidget *parent = nullptr)
            :
            QWidget(parent),
            moleculeWidget_(new MoleculeWidget(this)),
            idSpinBox_(new QSpinBox(this)),
            spinConnectionsCheckBox_(new QCheckBox("Spin Connections", this)),
            spinCorrelationsCheckBox_(new QCheckBox("Spin Correlations", this)),
            spinCorrelationSlider_(new QSlider(Qt::Orientation::Horizontal, this)),
            spinCorrelationSliderLabel_(new QLabel(this)) {
        auto splashScreen = createSplashScreen();

        loadData();

        auto hbox = new QHBoxLayout(this);
        auto gbox = new QGroupBox("Settings:");
        auto vbox = new QVBoxLayout(gbox);

        setLayout(hbox);
        resize(800, 800);
        hbox->addWidget(moleculeWidget_, Qt::AlignLeft);
        hbox->addWidget(gbox);
        gbox->setLayout(vbox);
        vbox->addWidget(idSpinBox_);
        vbox->addWidget(spinConnectionsCheckBox_);
        vbox->addWidget(spinCorrelationsCheckBox_);
        vbox->addWidget(spinCorrelationSlider_);
        vbox->addWidget(spinCorrelationSliderLabel_);


        QObject::connect(idSpinBox_, SIGNAL(valueChanged(int)),
                         this, SLOT(updateMoleculeWidget()));
        QObject::connect(spinConnectionsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(updateMoleculeWidget()));
        QObject::connect(spinCorrelationsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(updateMoleculeWidget()));
        QObject::connect(spinCorrelationSlider_, SIGNAL(valueChanged(int)),
                         this, SLOT(updateMoleculeWidget()));

        idSpinBox_->setRange(0, int(clusterCollection_.size() - 1));
        idSpinBox_->setSingleStep(1);
        idSpinBox_->setValue(0);

        spinCorrelationSlider_->setRange(0, 255);
        spinCorrelationSlider_->setSingleStep(1);
        spinCorrelationSlider_->setValue(191);

        spinCorrelationSliderLabel_->setFixedHeight(14);

        QTimer::singleShot(1000, splashScreen, SLOT(close()));
        QTimer::singleShot(1000, this, SLOT(show()));

        update();
    }

public slots:

    void updateMoleculeWidget() {

        auto spinCorrelationThreshold = double(spinCorrelationSlider_->value())/double(255);
        moleculeWidget_->setMolecule(
                atomsVector_,
                clusterCollection_[idSpinBox_->value()],
                spinConnectionsCheckBox_->checkState() == Qt::CheckState::Checked,
                spinCorrelationsCheckBox_->checkState() == Qt::CheckState::Checked,
                spinCorrelationThreshold
        );

        spinCorrelationSliderLabel_->setText(QString::number(spinCorrelationThreshold));
    };


private:
    MoleculeWidget *moleculeWidget_;
    QSpinBox *idSpinBox_;
    QCheckBox *spinConnectionsCheckBox_, *spinCorrelationsCheckBox_;
    QSlider *spinCorrelationSlider_;
    QLabel *spinCorrelationSliderLabel_;

    AtomsVector atomsVector_;
    std::vector<std::pair<std::vector<ElectronsVector>, YAML::Node>> clusterCollection_;


    QSplashScreen *createSplashScreen() {
        auto splashScreen = new QSplashScreen();

        QImage file(":inPsights.png");
        QPixmap qPixmap;
        qPixmap.convertFromImage(file);
        splashScreen->setPixmap(qPixmap);
        splashScreen->show();
        splashScreen->showMessage("Version 1.0.0", Qt::AlignRight, Qt::lightGray);

        return splashScreen;
    }

    void loadData(){

        auto dialog = new QFileDialog(this);
        dialog->setWindowTitle("Open results file");
        dialog->setFileMode(QFileDialog::FileMode::ExistingFile);
        dialog->setViewMode(QFileDialog::ViewMode::Detail);

        YAML::Node doc = YAML::LoadFile(dialog->getOpenFileName().toStdString());

        for(YAML::const_iterator it = doc["Clusters"].begin(); it != doc["Clusters"].end();++it) {
            clusterCollection_.emplace_back(
                    (*it)["Structures"].as<std::vector<ElectronsVector>>(),
                    (*it)["SpinCorrelations"]);
        }

        atomsVector_ = doc["Atoms"].as<AtomsVector>();

    }
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
