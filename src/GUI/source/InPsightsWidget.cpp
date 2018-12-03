//
// Created by heuer on 03.12.18.
//

#include <InPsightsWidget.h>
#include <QFileDialog>
#include <QString>

#include <QHBoxLayout>
#include <QGroupBox>
#include <QSpinBox>
#include <QSplashScreen>
#include <QTimer>

InPsightsWidget::InPsightsWidget(QWidget *parent)
            :
            QWidget(parent),
            moleculeWidget_(new MoleculeWidget(this)),
            bondsCheckBox_(new QCheckBox("Bonds", this)),
            spinConnectionsCheckBox_(new QCheckBox("Spin Connections", this)),
            spinCorrelationsCheckBox_(new QCheckBox("Spin Correlations", this)),
            spinCorrelationSlider_(new QSlider(Qt::Orientation::Horizontal, this)),
            spinCorrelationSliderLabel_(new QLabel(this)),
            maximaList_(new QListWidget(this))
    {
        auto splashScreen = createSplashScreen();

        loadData();

        auto gbox = new QGroupBox("Settings:");
        auto hbox = new QHBoxLayout(this);
        auto vboxOuter = new QVBoxLayout();
        auto vboxInner = new QVBoxLayout();


        setLayout(hbox);
        resize(1024, 768);
        hbox->addWidget(moleculeWidget_, Qt::AlignLeft);
        hbox->addLayout(vboxOuter);

        vboxOuter->addWidget(maximaList_);
        vboxOuter->addWidget(gbox);
        gbox->setLayout(vboxInner);

        vboxInner->addWidget(bondsCheckBox_);
        vboxInner->addWidget(spinConnectionsCheckBox_);
        vboxInner->addWidget(spinCorrelationsCheckBox_);
        vboxInner->addWidget(spinCorrelationSlider_);
        vboxInner->addWidget(spinCorrelationSliderLabel_);


        QObject::connect(maximaList_, SIGNAL(itemChanged(QListWidgetItem*)),
                         this, SLOT(selectedStructure(QListWidgetItem*)));

        QObject::connect(bondsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(onBondsChecked(int)));

        QObject::connect(spinConnectionsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(onSpinConnectionsChecked(int)));

        /*QObject::connect(maximaList_, SIGNAL(itemChanged(QListWidgetItem*)),
                this, SLOT(updateMoleculeWidget(QListWidgetItem*)));

        QObject::connect(spinCorrelationsCheckBox_, SIGNAL(stateChanged(int)),
                this, SLOT(updateMoleculeWidget()));
        QObject::connect(spinCorrelationSlider_, SIGNAL(valueChanged(int)),
                this, SLOT(updateMoleculeWidget()));*/

        spinCorrelationSlider_->setRange(0, 255);
        spinCorrelationSlider_->setSingleStep(1);
        spinCorrelationSlider_->setValue(191);

        spinCorrelationSliderLabel_->setFixedHeight(14);

        QTimer::singleShot(1000, splashScreen, SLOT(close()));
        QTimer::singleShot(1000, this, SLOT(show()));

        update();
        initialView();
    }

    void InPsightsWidget::selectedStructure(QListWidgetItem* item) {
        auto id = item->data(Qt::ItemDataRole::UserRole).toInt();
        auto data = clusterCollection_[id];

        if(item->checkState() == Qt::CheckState::Checked)
            moleculeWidget_->addElectronsVector(data.representativeStructure(), id);
        else
            moleculeWidget_->removeElectronsVector(id);
    };

    void InPsightsWidget::onBondsChecked(int stateId){
        moleculeWidget_->drawBonds(Qt::CheckState(stateId) == Qt::CheckState::Checked);
    }

    void InPsightsWidget::onSpinConnectionsChecked(int stateId){
        moleculeWidget_->drawSpinConnections(Qt::CheckState(stateId) == Qt::CheckState::Checked);
    }

    QSplashScreen *InPsightsWidget::createSplashScreen() {
        auto splashScreen = new QSplashScreen();

        QImage file(":inPsights.png");
        QPixmap qPixmap;
        qPixmap.convertFromImage(file);
        splashScreen->setPixmap(qPixmap);
        splashScreen->show();
        splashScreen->showMessage("Version 1.0.0", Qt::AlignRight, Qt::lightGray);

        return splashScreen;
    }

    void InPsightsWidget::loadData() {

        auto dialog = new QFileDialog(this);
        dialog->setWindowTitle("Open results file");
        dialog->setFileMode(QFileDialog::FileMode::ExistingFile);
        dialog->setViewMode(QFileDialog::ViewMode::Detail);

        YAML::Node doc = YAML::LoadFile(QFileDialog::getOpenFileName().toStdString());

        auto Vnn = doc["Vnn"];

        int id = 0;
        for(YAML::const_iterator it = doc["Clusters"].begin(); it != doc["Clusters"].end();++it) {
            //OneParticleEnergies::oneAtomEnergies(*it, Vnn);
            //OneParticleEnergies::oneElectronEnergies(*it);

            ClusterData clusterData = (*it).as<ClusterData>();

            auto text = QStringLiteral("Cluster %1 (%2)").arg(QString::number(id),QString::number(clusterData.N_));

            clusterCollection_.emplace_back(clusterData);

            auto item = new QListWidgetItem(text);
            item->setCheckState(Qt::CheckState::Unchecked);
            item->setData(Qt::ItemDataRole::UserRole, QVariant(id));

            maximaList_->addItem(item);
            id++;
        }
        moleculeWidget_->setSharedAtomsVector(doc["Atoms"].as<AtomsVector>());
    }

    void InPsightsWidget::initialView(){
        moleculeWidget_->drawAtoms();
        bondsCheckBox_->setCheckState(Qt::CheckState::Checked);
        maximaList_->item(0)->setCheckState(Qt::CheckState::Checked);
        spinConnectionsCheckBox_->setCheckState(Qt::CheckState::Checked);
    }


/*
void updateMoleculeWidget(QListWidgetItem* item) {
    auto spinCorrelationThreshold = double(spinCorrelationSlider_->value())/double(255);
    moleculeWidget_->setMolecule(
            atomsVector_,
            clusterCollection_[item->data(Qt::ItemDataRole::UserRole).toInt()],
            spinConnectionsCheckBox_->checkState() == Qt::CheckState::Checked,
            spinCorrelationsCheckBox_->checkState() == Qt::CheckState::Checked,
            spinCorrelationThreshold
    );
};

void updateMoleculeWidget() {

    auto spinCorrelationThreshold = double(spinCorrelationSlider_->value())/double(255);
    moleculeWidget_->setMolecule(
            atomsVector_,
            clusterCollection_[0],
            spinConnectionsCheckBox_->checkState() == Qt::CheckState::Checked,
            spinCorrelationsCheckBox_->checkState() == Qt::CheckState::Checked,
            spinCorrelationThreshold
    );

    spinCorrelationSliderLabel_->setText(QString::number(spinCorrelationThreshold));
};
*/