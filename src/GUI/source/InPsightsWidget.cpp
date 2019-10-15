/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <InPsightsWidget.h>
#include <QFileDialog>
#include <QString>

#include <QGridLayout>
#include <QGroupBox>
#include <QSplashScreen>
#include <QTimer>
#include <QHeaderView>
#include "IntegerSortedTreeWidgetItem.h"
#include <iterator>
#include <vector>
#include <ParticlesVector.h>
#include <SurfaceDataGenerator.h>
#include <spdlog/spdlog.h>


InPsightsWidget::InPsightsWidget(QWidget *parent, const std::string& filename)
        :
        QWidget(parent),
        filename_(filename),
        moleculeWidget(new MoleculeWidget(this)),
        maximaProcessingWidget(new MaximaProcessingWidget(this)), // TODO refator, should it be an additional window?
        atomsCheckBox(new QCheckBox("Atoms", this)),
        bondsCheckBox(new QCheckBox("Bonds", this)),
        axesCheckBox(new QCheckBox("Axes", this)),
        spinConnectionsCheckBox(new QCheckBox("Spin Connections", this)),
        spinCorrelationsCheckBox(new QCheckBox("Spin Correlations", this)),
        sedsCheckBox(new QCheckBox("SEDs", this)),
        maximaHullsCheckBox(new QCheckBox("Maxima Hulls", this)),
        spinCorrelationBox(new QDoubleSpinBox(this)),
        sedPercentageBox(new QDoubleSpinBox(this)),
        maximaList(new QTreeWidget(this)) {

    loadData();
    showSplashScreen();
    createWidget();
    connectSignals();
    initialView();
    show();
}

void InPsightsWidget::createWidget() {
    setWindowTitle("inPsights - Chemical insights from |Ψ|².");

    auto hbox = new QHBoxLayout(this);
    auto vboxOuter = new QVBoxLayout();
    auto vboxInner = new QVBoxLayout();
    auto gbox = new QGroupBox("Settings:");

    setLayout(hbox);

    resize(1280, 800);
    hbox->addWidget(moleculeWidget, 2);
    hbox->addLayout(vboxOuter, 1);

    // put into MaximaTreeWidget class
    auto headerLabels = QList<QString>({"ID", "Weight", "min(-ln|Ψ|²)", "max(-ln|Ψ|²)"});
    maximaList->setColumnCount(headerLabels.size());
    maximaList->setHeaderLabels(headerLabels);
    maximaList->header()->setStretchLastSection(false);
    maximaList->header()->setSectionResizeMode(0,QHeaderView::Stretch);
    maximaList->header()->setSectionResizeMode(1,QHeaderView::ResizeToContents);
    maximaList->header()->setSectionResizeMode(2,QHeaderView::ResizeToContents);
    maximaList->header()->setSectionResizeMode(3,QHeaderView::ResizeToContents);

    vboxOuter->addWidget(maximaList, 1);
    vboxOuter->addWidget(maximaProcessingWidget,1);
    vboxOuter->addWidget(gbox);
    gbox->setLayout(vboxInner);

    maximaList->setSortingEnabled(true);

    auto checkboxGrid = new QGridLayout();
    vboxInner->addLayout(checkboxGrid,1);
    checkboxGrid->addWidget(atomsCheckBox,0,0);
    checkboxGrid->addWidget(bondsCheckBox,1,0);
    checkboxGrid->addWidget(axesCheckBox,2,0);

    checkboxGrid->addWidget(spinConnectionsCheckBox,0,1);
    checkboxGrid->addWidget(maximaHullsCheckBox,1,1);
    checkboxGrid->addWidget(spinCorrelationsCheckBox,2,1);
    checkboxGrid->addWidget(spinCorrelationBox,2,2);
    checkboxGrid->addWidget(sedsCheckBox,3,1);
    checkboxGrid->addWidget(sedPercentageBox,3,2);

    setupSpinBoxes();
}

void InPsightsWidget::connectSignals() {
    connect(maximaList, &QTreeWidget::itemChanged,
            this, &InPsightsWidget::selectedStructure);

    connect(atomsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onAtomsChecked);

    connect(bondsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onBondsChecked);

    connect(axesCheckBox, &QCheckBox::stateChanged,
           this, &InPsightsWidget::onAxesChecked);

    connect(spinConnectionsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSpinConnectionsChecked);

    connect(spinCorrelationsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSpinCorrelationsChecked);

    connect(spinCorrelationsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSpinCorrelationsChecked);

    connect(spinCorrelationBox, qOverload<double>(&QDoubleSpinBox::valueChanged),
            this, &InPsightsWidget::onSpinCorrelationsBoxChanged);

    connect(maximaProcessingWidget, &MaximaProcessingWidget::atomsChecked,
            moleculeWidget, &MoleculeWidget::onAtomsChecked);
    connect(maximaProcessingWidget, &MaximaProcessingWidget::electronsChecked,
            moleculeWidget, &MoleculeWidget::onElectronsChecked);

    connect(maximaProcessingWidget, &MaximaProcessingWidget::atomsHighlighted,
            moleculeWidget, &MoleculeWidget::onAtomsHighlighted);
    connect(maximaProcessingWidget, &MaximaProcessingWidget::electronsHighlighted,
            moleculeWidget, &MoleculeWidget::onElectronsHighlighted);
}

void InPsightsWidget::setupSpinBoxes() {
    spinCorrelationBox->setRange(0.0,1.0);
    spinCorrelationBox->setSingleStep(0.01);
    spinCorrelationBox->setValue(1.0);

    sedPercentageBox->setRange(0.0,1.0);
    sedPercentageBox->setSingleStep(0.01);
    sedPercentageBox->setValue(0.2);
}

void InPsightsWidget::selectedStructure(QTreeWidgetItem *item, int column) {
    //auto maximaTreeWidgetItem = dynamic_cast<IntegerSortedTreeWidgetItem*>(item);

    if (column != 0)
        spdlog::critical("Column 0 expected but got {} ", column);

    auto id = item->data(0, Qt::ItemDataRole::UserRole).toList();
    auto clusterId = id[0].toInt();
    auto structureId = id[1].toInt();

    auto createQ = item->checkState(0) == Qt::CheckState::Checked;
    spdlog::info("Selected structure {0} from cluster {1} for {2}.", structureId, clusterId,
                  createQ ? "creation" : "deletion");

    if (createQ) {
        moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId], clusterId, structureId);
        maximaProcessingWidget->updateData(clusterCollection_[clusterId]);

        if(sedsCheckBox->checkState() == Qt::CheckState::Checked
        && moleculeWidget->activeSedsMap_.find(clusterId) == moleculeWidget->activeSedsMap_.end()) {
            if (clusterCollection_[clusterId].voxelCubes_.empty())
                spdlog::warn("Voxel cubes were not calculated.");
            else
                moleculeWidget->addSeds(clusterId, clusterCollection_, sedPercentageBox->value());
        }

        if(maximaHullsCheckBox->checkState() == Qt::CheckState::Checked
           && moleculeWidget->activeMaximaHullsMap_.find(clusterId) == moleculeWidget->activeMaximaHullsMap_.end()) {
            if (clusterCollection_[clusterId].exemplaricStructures_.empty())
                spdlog::warn("Voxel cubes were not calculated.");
            else
                moleculeWidget->addMaximaHulls(clusterId, clusterCollection_);
        }

    } else {
        moleculeWidget->removeElectronsVector(clusterId, structureId);

        if(moleculeWidget->activeSedsMap_.find(clusterId) != moleculeWidget->activeSedsMap_.end())
            moleculeWidget->removeSeds(clusterId);

        if(moleculeWidget->activeMaximaHullsMap_.find(clusterId) != moleculeWidget->activeMaximaHullsMap_.end())
            moleculeWidget->removeMaximaHulls(clusterId);
    }
    redrawSpinDecorations();
};

void InPsightsWidget::redrawSpinDecorations() {
    if(spinConnectionsCheckBox->checkState() == Qt::CheckState::Checked){
        onSpinConnectionsChecked(Qt::Unchecked);
        onSpinConnectionsChecked(Qt::Checked);
    }
    if(spinCorrelationsCheckBox->checkState() == Qt::Checked){
        onSpinCorrelationsChecked(Qt::Unchecked);
        onSpinCorrelationsChecked(Qt::Checked);
    }
}

void InPsightsWidget::onAtomsChecked(int stateId) {
    moleculeWidget->drawAtoms(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onBondsChecked(int stateId) {
    moleculeWidget->drawBonds(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onAxesChecked(int stateId) {
    moleculeWidget->drawAxes(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onSpinConnectionsChecked(int stateId) {
    moleculeWidget->drawSpinConnections(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onSpinCorrelationsChecked(int stateId) {
    moleculeWidget->drawSpinCorrelations(Qt::CheckState(stateId) == Qt::CheckState::Checked,
                                         clusterCollection_, spinCorrelationBox->value());
}

void InPsightsWidget::onSpinCorrelationsBoxChanged(double value) {
    if (spinCorrelationsCheckBox->checkState() == Qt::CheckState::Checked) {
        onSpinCorrelationsChecked(Qt::CheckState::Unchecked); //TODO ugly, create update() function in SpinCorrelation3D and make it accessible
        onSpinCorrelationsChecked(Qt::CheckState::Checked);
    }
}

#include <QSize>
void InPsightsWidget::showSplashScreen() {
    auto splashScreen = new QSplashScreen();
    auto pixmap = QPixmap(":inPsights.png").scaledToWidth(400, Qt::TransformationMode::SmoothTransformation);

    splashScreen->setPixmap(pixmap);
    splashScreen->show();
    splashScreen->showMessage(
            "Version 0.1.0 (pre-release)\n"
            "Copyright © 2016-2019  Michael A. Heuer.", Qt::AlignBottom, Qt::gray);

    splashScreen->finish(this);
}

void InPsightsWidget::loadData() {

    if(filename_.empty()) {
        filename_ = QFileDialog::getOpenFileName(this,
                                                QString("Open results file"),
                                                QDir::currentPath(),
                                                QString("YAML files (*.yml *.yaml *.json)")).toStdString();
    }

    moleculeWidget->infoText_->setText(filename_.c_str());

    YAML::Node doc = YAML::LoadAllFromFile(filename_)[1]; // load results
    auto atoms = doc["Atoms"].as<AtomsVector>();

    auto nElectrons = doc["Clusters"][0]["Structures"][0].as<ElectronsVector>().numberOfEntities();

    auto EnStats = doc["En"].as<VectorStatistics>();
    maximaProcessingWidget->setAtomEnergies(EnStats);
    maximaProcessingWidget->setAtomsVector(atoms);
    maximaProcessingWidget->initializeTreeItems(maximaProcessingWidget->atomsTreeWidget(), int(atoms.numberOfEntities()));
    maximaProcessingWidget->initializeTreeItems(maximaProcessingWidget->electronsTreeWidget(), int(nElectrons));


    moleculeWidget->setSharedAtomsVector(atoms);


    for (int clusterId = 0; clusterId < static_cast<int>(doc["Clusters"].size()); ++clusterId) {

        ClusterData clusterData = doc["Clusters"][clusterId].as<ClusterData>();

        clusterCollection_.emplace_back(clusterData);
        auto item = new IntegerSortedTreeWidgetItem(
                maximaList, {QString::number(clusterId),
                 QString::number(1.0 * clusterData.N_ / doc["NSamples"].as<unsigned>(), 'f', 4),
                 QString::number(clusterData.valueStats_.cwiseMin()[0], 'f', 3),
                 QString::number(clusterData.valueStats_.cwiseMax()[0], 'f', 3)});

        item->setCheckState(0, Qt::CheckState::Unchecked);

        QList<QVariant> id = {clusterId, 0};
        item->setData(0, Qt::ItemDataRole::UserRole, id);

        auto structures = doc["Clusters"][clusterId]["Structures"];

        for (int structureId = 1; structureId < static_cast<int>(structures.size()); ++structureId) {
            auto subItem = new IntegerSortedTreeWidgetItem(item, QStringList({QString::number(structureId)}));
            subItem->setCheckState(0, Qt::CheckState::Unchecked);

            id = {clusterId, structureId};
            subItem->setData(0, Qt::ItemDataRole::UserRole, id);
            item->addChild(subItem);
        }

        maximaList->addTopLevelItem(item);
        for (int i = 0; i < maximaList->columnCount(); ++i) {
            maximaList->resizeColumnToContents(i);
        }
    }

}

void InPsightsWidget::initialView() {
    maximaList->sortItems(0,Qt::SortOrder::AscendingOrder);
    atomsCheckBox->setCheckState(Qt::CheckState::Checked);
    bondsCheckBox->setCheckState(Qt::CheckState::Checked);
    axesCheckBox->setCheckState(Qt::CheckState::Unchecked);
    maximaList->topLevelItem(0)->setCheckState(0, Qt::CheckState::Checked);
}
