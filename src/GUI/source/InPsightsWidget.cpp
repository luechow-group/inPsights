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
#include <QTreeWidgetItem>
#include <iterator>
#include <vector>

InPsightsWidget::InPsightsWidget(QWidget *parent)
            :
            QWidget(parent),
            moleculeWidget_(new MoleculeWidget(this)),
            atomsCheckBox_(new QCheckBox("Atoms", this)),
            bondsCheckBox_(new QCheckBox("Bonds", this)),
            spinConnectionsCheckBox_(new QCheckBox("Spin Connections", this)),
            spinCorrelationsCheckBox_(new QCheckBox("Spin Correlations", this)),
            spinCorrelationSlider_(new QSlider(Qt::Orientation::Horizontal, this)),
            spinCorrelationSliderLabel_(new QLabel(this)),
            maximaList_(new QTreeWidget(this))
    {
        auto splashScreen = createSplashScreen();
        setWindowIcon(QIcon(":inPsightsIcon.png"));
        setWindowTitle("inPsights - Chemical insights from |Ψ|².");

        loadData();

        auto gbox = new QGroupBox("Settings:");
        auto hbox = new QHBoxLayout(this);
        auto vboxOuter = new QVBoxLayout();
        auto vboxInner = new QVBoxLayout();

        setLayout(hbox);
        resize(1024, 768);
        hbox->addWidget(moleculeWidget_, Qt::AlignLeft);
        hbox->addLayout(vboxOuter);

        auto headerLabels = QList<QString>({"ID","N", "min(-ln(|Ψ|²))","max(-ln(|Ψ|²))"});
        maximaList_->setColumnCount(headerLabels.size());
        maximaList_->setHeaderLabels(headerLabels);
        vboxOuter->addWidget(maximaList_);
        vboxOuter->addWidget(gbox);
        gbox->setLayout(vboxInner);

        maximaList_->setFixedWidth(300);
        vboxInner->addWidget(atomsCheckBox_);
        vboxInner->addWidget(bondsCheckBox_);
        vboxInner->addWidget(spinConnectionsCheckBox_);
        vboxInner->addWidget(spinCorrelationsCheckBox_);
        vboxInner->addWidget(spinCorrelationSlider_);
        vboxInner->addWidget(spinCorrelationSliderLabel_);


        QObject::connect(maximaList_, SIGNAL(itemChanged(QTreeWidgetItem*,int)),
                         this, SLOT(selectedStructure(QTreeWidgetItem*, int)));

        QObject::connect(atomsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(onAtomsChecked(int)));

        QObject::connect(bondsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(onBondsChecked(int)));

        QObject::connect(spinConnectionsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(onSpinConnectionsChecked(int)));

        QObject::connect(spinCorrelationsCheckBox_, SIGNAL(stateChanged(int)),
                         this, SLOT(onSpinCorrelationsChecked(int)));


        //QObject::connect(spinCorrelationSlider_, SIGNAL(valueChanged(int)),
        //                 this, SLOT(updateMoleculeWidget()));

        spinCorrelationSlider_->setRange(0, 255);
        spinCorrelationSlider_->setSingleStep(1);
        spinCorrelationSlider_->setValue(191);

        spinCorrelationSliderLabel_->setFixedHeight(14);

        QTimer::singleShot(1000, splashScreen, SLOT(close()));
        QTimer::singleShot(1000, this, SLOT(show()));

        update();
        initialView();
    }

    void InPsightsWidget::selectedStructure(QTreeWidgetItem* item, int column) {

        if(column == 0) {
            auto clusterId = item->data(0, Qt::ItemDataRole::UserRole).toInt();
            auto structureId = item->data(1, Qt::ItemDataRole::UserRole).toInt();

            if (item->checkState(0) == Qt::CheckState::Checked)
                moleculeWidget_->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId],
                        clusterId, structureId);
            else
                moleculeWidget_->removeElectronsVector(clusterId,structureId);
        } else {
            std::cout << "wrong column" << std::endl; // TODO treat properly
        }
    };

    void InPsightsWidget::onAtomsChecked(int stateId){
        moleculeWidget_->drawAtoms(Qt::CheckState(stateId) == Qt::CheckState::Checked);
    }

    void InPsightsWidget::onBondsChecked(int stateId){
        moleculeWidget_->drawBonds(Qt::CheckState(stateId) == Qt::CheckState::Checked);
    }

    void InPsightsWidget::onSpinConnectionsChecked(int stateId){
        moleculeWidget_->drawSpinConnections(Qt::CheckState(stateId) == Qt::CheckState::Checked);
    }

    void InPsightsWidget::onSpinCorrelationsChecked(int stateId) {
        std::cout << "WRONG STATS ARE USED" << std::endl;
        moleculeWidget_->drawSpinCorrelations(Qt::CheckState(stateId) == Qt::CheckState::Checked,
                clusterCollection_[0].SeeStats_, //TODO !!! WRONG ONE
                spinCorrelationSlider_->value()
                );
    }

    void InPsightsWidget::onSpinCorrelationsSliderChanged(int value) {
        //if(spinConnectionsCheckBox_->checkState() == Qt::CheckState::Checked)
        //    //TODO CONTINUE HERE

    }

QSplashScreen *InPsightsWidget::createSplashScreen() {
        auto splashScreen = new QSplashScreen();

        splashScreen->setPixmap(QPixmap(":inPsights.png"));
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

        for (int clusterId = 0; clusterId < static_cast<int>(doc["Clusters"].size()); ++clusterId) {
            //OneParticleEnergies::oneAtomEnergies(*it, Vnn);
            //OneParticleEnergies::oneElectronEnergies(*it);

            ClusterData clusterData = doc["Clusters"][clusterId].as<ClusterData>();

            clusterCollection_.emplace_back(clusterData);
            auto item = new QTreeWidgetItem(QList<QString>({
                QString::number(clusterId),
                QString::number(clusterData.N_),
                QString::number(clusterData.valueStats_.cwiseMin()[0]),
                QString::number(clusterData.valueStats_.cwiseMax()[0])
            }));

            item->setCheckState(0,Qt::CheckState::Unchecked);
            item->setData(0, Qt::ItemDataRole::UserRole, QVariant(clusterId));
            item->setData(1, Qt::ItemDataRole::UserRole, QVariant(0)); // representative structure

            auto structures = doc["Clusters"][clusterId]["Structures"];

            for (int structure = 1; structure < static_cast<int>(structures.size()); ++structure) {
                auto subItem = new QTreeWidgetItem(QList<QString>({QString::number(structure)}));
                subItem->setCheckState(0,Qt::CheckState::Unchecked);
                subItem->setData(0, Qt::ItemDataRole::UserRole, QVariant(clusterId));
                subItem->setData(1, Qt::ItemDataRole::UserRole, QVariant(structure));
                item->addChild(subItem);
            }

            maximaList_->addTopLevelItem(item);
        }
        moleculeWidget_->setSharedAtomsVector(doc["Atoms"].as<AtomsVector>());
    }

    void InPsightsWidget::initialView(){
        maximaList_->resizeColumnToContents(0);
        maximaList_->resizeColumnToContents(1);
        maximaList_->resizeColumnToContents(2);
        maximaList_->resizeColumnToContents(3);
        atomsCheckBox_->setCheckState(Qt::CheckState::Checked);
        bondsCheckBox_->setCheckState(Qt::CheckState::Checked);
        maximaList_->topLevelItem(0)->setCheckState(0,Qt::CheckState::Checked);
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