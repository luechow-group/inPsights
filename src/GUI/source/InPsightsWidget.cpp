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
        console(spdlog::get(Logger::name)),
        moleculeWidget(new MoleculeWidget(this)),
        atomsCheckBox(new QCheckBox("Atoms", this)),
        bondsCheckBox(new QCheckBox("Bonds", this)),
        spinConnectionsCheckBox(new QCheckBox("Spin Connections", this)),
        spinCorrelationsCheckBox(new QCheckBox("Spin Correlations", this)),
        spinCorrelationSlider(new QSlider(Qt::Orientation::Horizontal, this)),
        spinCorrelationSliderLabel(new QLabel(this)),
        maximaList(new QTreeWidget(this)) {
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
    hbox->addWidget(moleculeWidget, Qt::AlignLeft);
    hbox->addLayout(vboxOuter);

    auto headerLabels = QList<QString>({"ID", "N", "min(-ln(|Ψ|²))", "max(-ln(|Ψ|²))"});
    maximaList->setColumnCount(headerLabels.size());
    maximaList->setHeaderLabels(headerLabels);
    vboxOuter->addWidget(maximaList);
    vboxOuter->addWidget(gbox);
    gbox->setLayout(vboxInner);

    maximaList->setFixedWidth(300);
    vboxInner->addWidget(atomsCheckBox);
    vboxInner->addWidget(bondsCheckBox);
    vboxInner->addWidget(spinConnectionsCheckBox);
    vboxInner->addWidget(spinCorrelationsCheckBox);
    vboxInner->addWidget(spinCorrelationSlider);
    vboxInner->addWidget(spinCorrelationSliderLabel);


    QObject::connect(maximaList, SIGNAL(itemChanged(QTreeWidgetItem * , int)),
                     this, SLOT(selectedStructure(QTreeWidgetItem * , int)));

    QObject::connect(atomsCheckBox, SIGNAL(stateChanged(int)),
                     this, SLOT(onAtomsChecked(int)));

    QObject::connect(bondsCheckBox, SIGNAL(stateChanged(int)),
                     this, SLOT(onBondsChecked(int)));

    QObject::connect(spinConnectionsCheckBox, SIGNAL(stateChanged(int)),
                     this, SLOT(onSpinConnectionsChecked(int)));

    QObject::connect(spinCorrelationsCheckBox, SIGNAL(stateChanged(int)),
                     this, SLOT(onSpinCorrelationsChecked(int)));


    QObject::connect(spinCorrelationSlider, SIGNAL(valueChanged(int)),
                     this, SLOT(onSpinCorrelationsSliderChanged(int)));

    spinCorrelationSlider->setRange(0, 255);
    spinCorrelationSlider->setSingleStep(1);
    spinCorrelationSlider->setValue(int(0.75*255));

    spinCorrelationSliderLabel->setFixedHeight(14);

    QTimer::singleShot(1000, splashScreen, SLOT(close()));
    QTimer::singleShot(1000, this, SLOT(show()));

    update();
    initialView();
}

void InPsightsWidget::selectedStructure(QTreeWidgetItem *item, int column) {

    if (column != 0)
        console->critical("Column 0 expected but got {} ", column);

    auto clusterId = item->data(0, Qt::ItemDataRole::UserRole).toInt();
    auto structureId = item->data(1, Qt::ItemDataRole::UserRole).toInt();

    auto createQ = item->checkState(0) == Qt::CheckState::Checked;
    console->info("Selected structure {1} from cluster {0} for {2}.", clusterId, structureId,
                  createQ ? "creation" : "deletion");

    if (createQ)
        moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId], clusterId,
                                           structureId);
    else
        moleculeWidget->removeElectronsVector(clusterId, structureId);
};

void InPsightsWidget::onAtomsChecked(int stateId) {
    moleculeWidget->drawAtoms(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onBondsChecked(int stateId) {
    moleculeWidget->drawBonds(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onSpinConnectionsChecked(int stateId) {
    moleculeWidget->drawSpinConnections(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onSpinCorrelationsChecked(int stateId) {
    console->critical("ATTENTION: Currently, the same correlation data is plotted for all structures.");
    moleculeWidget->drawSpinCorrelations(Qt::CheckState(stateId) == Qt::CheckState::Checked,
                                         clusterCollection_[0].SeeStats_, //TODO !!! WRONG ONE
                                         spinCorrelationSlider->value()
    );
}

void InPsightsWidget::onSpinCorrelationsSliderChanged(int value) {
    if(spinCorrelationsCheckBox->checkState() == Qt::CheckState::Checked){
        onSpinCorrelationsChecked(Qt::CheckState::Unchecked); //TODO ugly, create update() function in SpinCorrelation3D and make it accessible
        onSpinCorrelationsChecked(Qt::CheckState::Checked);
    }
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
        auto item = new QTreeWidgetItem(QList<QString>({QString::number(clusterId),
                                                        QString::number(clusterData.N_),
                                                        QString::number(clusterData.valueStats_.cwiseMin()[0]),
                                                        QString::number(clusterData.valueStats_.cwiseMax()[0])}));

        item->setCheckState(0, Qt::CheckState::Unchecked);
        item->setData(0, Qt::ItemDataRole::UserRole, QVariant(clusterId));
        item->setData(1, Qt::ItemDataRole::UserRole, QVariant(0)); // representative structure

        auto structures = doc["Clusters"][clusterId]["Structures"];

        for (int structure = 1; structure < static_cast<int>(structures.size()); ++structure) {
            auto subItem = new QTreeWidgetItem(QList<QString>({QString::number(structure)}));
            subItem->setCheckState(0, Qt::CheckState::Unchecked);
            subItem->setData(0, Qt::ItemDataRole::UserRole, QVariant(clusterId));
            subItem->setData(1, Qt::ItemDataRole::UserRole, QVariant(structure));
            item->addChild(subItem);
        }

        maximaList->addTopLevelItem(item);
    }
    moleculeWidget->setSharedAtomsVector(doc["Atoms"].as<AtomsVector>());
}

void InPsightsWidget::initialView() {
    maximaList->resizeColumnToContents(0);
    maximaList->resizeColumnToContents(1);
    maximaList->resizeColumnToContents(2);
    maximaList->resizeColumnToContents(3);
    atomsCheckBox->setCheckState(Qt::CheckState::Checked);
    bondsCheckBox->setCheckState(Qt::CheckState::Checked);
    maximaList->topLevelItem(0)->setCheckState(0, Qt::CheckState::Checked);
    spinConnectionsCheckBox->setCheckState(Qt::CheckState::Checked);
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