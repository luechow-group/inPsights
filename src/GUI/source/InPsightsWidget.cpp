// Copyright (C) 2018-2021 Michael Heuer.
// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <InPsightsWidget.h>
#include <QFileDialog>
#include <QString>
#include <QGridLayout>
#include <QGroupBox>
#include <QSplashScreen>
#include <QTimer>
#include <QHeaderView>
#include "IntegerSortedTreeWidgetItem.h"
#include <BuildInfo.h>
#include <vector>
#include <ParticlesVector.h>
#include "CameraSettings.h"
#include <spdlog/spdlog.h>


InPsightsWidget::InPsightsWidget(QWidget *parent, const std::string& filename)
        :
        QWidget(parent),
        filename_(filename),
        moleculeWidget(new MoleculeWidget(this)),
        //maximaProcessingWidget(new MaximaProcessingWidget(this)),
        atomsCheckBox(new QCheckBox("Nuclei", this)),
        bondsCheckBox(new QCheckBox("Bonds", this)),
        axesCheckBox(new QCheckBox("Axes", this)),
        sampleAverageCheckBox(new QCheckBox("Sample Average", this)),
        spinCorrelationsCheckBox(new QCheckBox("Spin Correlations", this)),
        sedsExportButton(new QPushButton("export SEDs", this)),
        sedsCheckBox(new QCheckBox("SEDs", this)),
        maximaHullsCheckBox(new QCheckBox("Maxima Hulls", this)),
        plotAllCheckBox(new QCheckBox("All of Cluster", this)),
        coloredCheckBox(new QCheckBox("Multicolored", this)),
        spinCorrelationBox(new QDoubleSpinBox(this)),
        sedPercentageBox(new QDoubleSpinBox(this)),
        atom1Box(new QSpinBox(this)),
        atom2Box(new QSpinBox(this)),
        electron1Box(new QSpinBox(this)),
        electron2Box(new QSpinBox(this)),
        maximaList(new QTreeWidget(this)),
        probabilitySum(new QLabel(this)),
        deselectAllButton(new QPushButton("deselect all", this))
        {

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
    auto vboxInner2 = new QVBoxLayout();
    auto selectorBox = new QGroupBox("Particle highlighting:");
    auto settingsBox = new QGroupBox("Settings:");

    setLayout(hbox);

    resize(1424, 801);
    hbox->addWidget(moleculeWidget, 5);
    hbox->addLayout(vboxOuter, 1);

    // put into MaximaTreeWidget class
    auto headerLabels = QList<QString>({"ID", "P", "min(Φ)", "max(Φ)"});
    maximaList->setColumnCount(headerLabels.size());
    maximaList->setHeaderLabels(headerLabels);
    maximaList->header()->setStretchLastSection(false);

    maximaList->headerItem()->setToolTip(0,QString("IDs sorted by potential Φ"));
    maximaList->headerItem()->setToolTip(1,QString("Probability"));
    maximaList->headerItem()->setToolTip(2,QString("Φ = -ħ/2m ln|Ψ|²"));
    maximaList->headerItem()->setToolTip(3,QString("Φ = -ħ/2m ln|Ψ|²"));

    maximaList->header()->setSectionResizeMode(0,QHeaderView::Stretch);
    maximaList->header()->setSectionResizeMode(1,QHeaderView::ResizeToContents);
    maximaList->header()->setSectionResizeMode(2,QHeaderView::ResizeToContents);
    maximaList->header()->setSectionResizeMode(3,QHeaderView::ResizeToContents);

    vboxOuter->addWidget(maximaList, 1);

    auto lineGrid = new QGridLayout();
    vboxOuter->addLayout(lineGrid,1);
    lineGrid->addWidget(probabilitySum,0,0);
    lineGrid->addWidget(deselectAllButton,0,1);

    vboxOuter->addWidget(selectorBox);
    selectorBox->setLayout(vboxInner2);

    auto selectorGrid = new QGridLayout();
    vboxInner2->addLayout(selectorGrid,1);
    selectorGrid->addWidget(new QLabel("Nuclei"),0,0);
    selectorGrid->addWidget(atom1Box,1,0);
    selectorGrid->addWidget(atom2Box,1,1);
    selectorGrid->addWidget(new QLabel("Electrons"),2,0);
    selectorGrid->addWidget(electron1Box,3,0);
    selectorGrid->addWidget(electron2Box,3,1);

    //vboxOuter->addWidget(maximaProcessingWidget,1);
    vboxOuter->addWidget(settingsBox);
    settingsBox->setLayout(vboxInner);

    maximaList->setSortingEnabled(true);

    auto checkboxGrid = new QGridLayout();
    vboxInner->addLayout(checkboxGrid,1);

    // first row
    checkboxGrid->addWidget(atomsCheckBox,0,0);
    checkboxGrid->addWidget(bondsCheckBox,0,1);

    //second row
    checkboxGrid->addWidget(axesCheckBox,1,0);
    checkboxGrid->addWidget(sampleAverageCheckBox,1,1);


    //third row
    checkboxGrid->addWidget(coloredCheckBox,2,0);
    checkboxGrid->addWidget(maximaHullsCheckBox,2,1);

    //fourth row
    checkboxGrid->addWidget(sedsExportButton,3,0);
    checkboxGrid->addWidget(plotAllCheckBox,3,1);

    //fifth row
    checkboxGrid->addWidget(sedsCheckBox,4,0);
    checkboxGrid->addWidget(sedPercentageBox,4,1);

    //sixth row
    checkboxGrid->addWidget(spinCorrelationsCheckBox,5,0);
    checkboxGrid->addWidget(spinCorrelationBox,5,1);

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

    connect(spinCorrelationsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSpinCorrelationsChecked);

    connect(spinCorrelationsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSpinCorrelationsChecked);

    connect(plotAllCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onPlotAllChecked);

    connect(spinCorrelationBox, qOverload<double>(&QDoubleSpinBox::valueChanged),
            this, &InPsightsWidget::onSpinCorrelationsBoxChanged);

    connect(atom1Box, qOverload<int>(&QSpinBox::valueChanged),
            this, &InPsightsWidget::onAtom1BoxChanged);

    connect(atom2Box, qOverload<int>(&QSpinBox::valueChanged),
            this, &InPsightsWidget::onAtom2BoxChanged);

    connect(electron1Box, qOverload<int>(&QSpinBox::valueChanged),
            this, &InPsightsWidget::onElectron1BoxChanged);

    connect(electron2Box, qOverload<int>(&QSpinBox::valueChanged),
            this, &InPsightsWidget::onElectron2BoxChanged);

    connect(sedsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSedChecked);

    connect(coloredCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::updateSelectedStructures);

    connect(maximaHullsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::updateSelectedStructures);

    connect(sampleAverageCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::updateSelectedStructures);

    /*
    connect(maximaProcessingWidget, &MaximaProcessingWidget::atomsChecked,
            moleculeWidget, &MoleculeWidget::onAtomsChecked);
    connect(maximaProcessingWidget, &MaximaProcessingWidget::electronsChecked,
            moleculeWidget, &MoleculeWidget::onElectronsChecked);

    connect(maximaProcessingWidget, &MaximaProcessingWidget::atomsHighlighted,
            moleculeWidget, &MoleculeWidget::onAtomsHighlighted);
    connect(maximaProcessingWidget, &MaximaProcessingWidget::electronsHighlighted,
            moleculeWidget, &MoleculeWidget::onElectronsHighlighted);
    */
    connect(sedsExportButton, &QPushButton::clicked, this, &InPsightsWidget::onSedsExport);
    connect(deselectAllButton, &QPushButton::clicked, this, &InPsightsWidget::onDeselectAll);
}

void InPsightsWidget::setupSpinBoxes() {
    spinCorrelationBox->setRange(0.0,1.0);
    spinCorrelationBox->setSingleStep(0.01);
    spinCorrelationBox->setValue(1.0);

    sedPercentageBox->setRange(0.0,1.0);
    sedPercentageBox->setSingleStep(0.01);
    sedPercentageBox->setValue(0.2);

    atom1Box->setRange(-1,moleculeWidget->getAtomsNumber()-1);
    atom1Box->setSingleStep(1);
    atom1Box->setValue(-1);

    atom2Box->setRange(-1,moleculeWidget->getAtomsNumber()-1);
    atom2Box->setSingleStep(1);
    atom2Box->setValue(-1);

    electron1Box->setRange(-1,clusterCollection_[0].representativeStructure().numberOfEntities()-1);
    electron1Box->setSingleStep(1);
    electron1Box->setValue(-1);

    electron2Box->setRange(-1,clusterCollection_[0].representativeStructure().numberOfEntities()-1);
    electron2Box->setSingleStep(1);
    electron2Box->setValue(-1);
}

void InPsightsWidget::selectedStructure(QTreeWidgetItem *item, int column) {
    //auto maximaTreeWidgetItem = dynamic_cast<IntegerSortedTreeWidgetItem*>(item);

    if (column != 0)
        spdlog::critical("Column 0 expected but got {} ", column);

    auto id = item->data(0, Qt::ItemDataRole::UserRole).toList();
    auto clusterId = id[0].toInt();
    auto secondId = id[1].toInt();
    auto structureId = secondId;
    if (structureId == -1){  // if secondId == -1, the top level checkbox is marked
        structureId = 0;  // taking structure 0 to visualize the top level cluster
    }

    auto createQ = item->checkState(0) == Qt::CheckState::Checked;

    if (createQ) {
        auto sampleAverage = clusterCollection_[clusterId].sampleAverage_;

        if (sampleAverageCheckBox->checkState() == Qt::CheckState::Checked &&
            sampleAverage.numberOfEntities() > 0) {
            moleculeWidget->addElectronsVector(sampleAverage, clusterId, secondId,
                                               coloredCheckBox->checkState() == Qt::Checked);
        } else if (sampleAverageCheckBox->checkState() == Qt::CheckState::Checked
                   && sampleAverage.numberOfEntities() == 0) {
            moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId],
                                               clusterId, secondId, coloredCheckBox->checkState() == Qt::Checked);
            spdlog::warn("Sample averaged vectors were not calculated. Plotting the first maximum instead...");
        } else {
            moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId],
                                               clusterId, secondId, coloredCheckBox->checkState() == Qt::Checked);
        }

        //maximaProcessingWidget->updateData(clusterCollection_[clusterId]);

        if (sedsCheckBox->checkState() == Qt::CheckState::Checked
            && moleculeWidget->activeSedsMap_.find(clusterId) == moleculeWidget->activeSedsMap_.end()) {
            moleculeWidget->addSeds(clusterId, structureId, clusterCollection_, sedPercentageBox->value());
        }

        if (maximaHullsCheckBox->checkState() == Qt::CheckState::Checked
            &&
            moleculeWidget->activeMaximaHullsMap_.find(clusterId) == moleculeWidget->activeMaximaHullsMap_.end()) {
            if (clusterCollection_[clusterId].exemplaricStructures_.empty())
                spdlog::warn("Voxel cubes were not calculated.");
            else
                moleculeWidget->addMaximaHulls(clusterId, clusterCollection_);
        }

        onElectron1BoxChanged(electron1Box->value());
        onElectron2BoxChanged(electron2Box->value());
    } else {
        moleculeWidget->removeElectronsVector(clusterId, secondId);

        if (moleculeWidget->activeSedsMap_.find(clusterId) != moleculeWidget->activeSedsMap_.end())
            moleculeWidget->removeSeds(clusterId);

        if (moleculeWidget->activeMaximaHullsMap_.find(clusterId) != moleculeWidget->activeMaximaHullsMap_.end())
            moleculeWidget->removeMaximaHulls(clusterId);
    }
    redrawSpinDecorations();
    probabilitySum->setText(QString("Σ=") + QString::number(sumProbabilities(), 'f', 4));

};

void InPsightsWidget::updateSelectedStructures(int) {
    auto* root = maximaList->invisibleRootItem();

    // iterate over topLevelItems
    for (int i = 0; i < root->childCount(); ++i) {
        if(root->child(i)->checkState(0) == Qt::Checked) {
            root->child(i)->setCheckState(0, Qt::Unchecked);
            root->child(i)->setCheckState(0, Qt::Checked);
        }
        // iterate over childs of topLevelItem i
        for (int j = 0; j<root->child(i)->childCount(); ++j) {
            if(root->child(i)->child(j)->checkState(0) == Qt::Checked) {
                root->child(i)->child(j)->setCheckState(0, Qt::Unchecked);
                root->child(i)->child(j)->setCheckState(0, Qt::Checked);
            }
        }
    }
};

void InPsightsWidget::onSedChecked(int stateId) {
    // guessing that all voxelCubes are empty, if first is empty. Should be fine
    if(clusterCollection_[0].voxelCubes_.empty()) {
        if (Qt::CheckState(stateId) == Qt::CheckState::Checked) {
            spdlog::warn("Voxel cubes were not calculated.");
            sedsCheckBox->setCheckState(Qt::CheckState::Unchecked);
        }
    }
    else {
        updateSelectedStructures(Qt::CheckState::Checked);
    }
}

void InPsightsWidget::onPlotAllChecked(int stateId)  {
    auto plotAllQ = Qt::CheckState(stateId);

    auto* root = maximaList->invisibleRootItem();

    // iterate over topLevelItems
    for (int i = 0; i < root->childCount(); ++i) {
        if(root->child(i)->checkState(0) == Qt::Checked) {

            // iterate over childs of topLevelItem i
            for (int j = 0; j < root->child(i)->childCount(); ++j)
                root->child(i)->child(j)->setCheckState(0, plotAllQ);
        }
    }
}

void InPsightsWidget::redrawSpinDecorations() {
    onSpinCorrelationsChecked(spinCorrelationsCheckBox->checkState());
}

void InPsightsWidget::onAtomsChecked(int stateId) {
    bondsCheckBox->setCheckState(Qt::CheckState::Unchecked);
    moleculeWidget->drawAtoms(Qt::CheckState(stateId) == Qt::CheckState::Checked);
    if (Qt::CheckState(stateId) == Qt::CheckState::Checked) {
        onAtom1BoxChanged(atom1Box->value());
        onAtom2BoxChanged(atom2Box->value());
    }
}

void InPsightsWidget::onBondsChecked(int stateId) {
    if (atomsCheckBox->isChecked()) {
        moleculeWidget->drawBonds(Qt::CheckState(stateId) == Qt::CheckState::Checked);
    }
    else {
        if (Qt::CheckState(stateId) == Qt::CheckState::Checked) {
            bondsCheckBox->setCheckState(Qt::CheckState::Unchecked);
        }
    }
}

void InPsightsWidget::onAxesChecked(int stateId) {
    moleculeWidget->drawAxes(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onSpinCorrelationsChecked(int stateId) {
    moleculeWidget->deleteSpinCorrelations();
    if(Qt::CheckState(stateId) == Qt::Checked)
        moleculeWidget->drawSpinCorrelations(clusterCollection_, spinCorrelationBox->value(),true);
    else if(Qt::CheckState(stateId) == Qt::PartiallyChecked)
        moleculeWidget->drawSpinCorrelations(clusterCollection_, spinCorrelationBox->value(),false);
}

void InPsightsWidget::onSpinCorrelationsBoxChanged(double value) {
    redrawSpinDecorations();
}

void InPsightsWidget::onAtom1BoxChanged(int value) {
    moleculeWidget->onAtomsHighlighted(value);
}

void InPsightsWidget::onAtom2BoxChanged(int value) {
    moleculeWidget->onAtomsChecked(value);
}

void InPsightsWidget::onElectron1BoxChanged(int value) {
    moleculeWidget->onElectronsHighlighted(value);
}

void InPsightsWidget::onElectron2BoxChanged(int value) {
    moleculeWidget->onElectronsChecked(value);
}

void InPsightsWidget::showSplashScreen() {
    auto splashScreen = new QSplashScreen();
    auto pixmap = QPixmap(":inPsights.png").scaledToWidth(450, Qt::TransformationMode::SmoothTransformation);

    splashScreen->setPixmap(pixmap);
    splashScreen->show();

    std::string message = inPsights::version() + "\n"\
                          "Copyright © 2016-2021  Michael A. Heuer.\n"\
                          "Copyright © 2018-2022  Leonard Reuter.";

    splashScreen->showMessage(message.c_str(), Qt::AlignBottom, Qt::gray);

    splashScreen->finish(this);
}

void InPsightsWidget::loadData() {

    if(filename_.empty()) {
        filename_ = QFileDialog::getOpenFileName(this,
                                                QString("Open results file"),
                                                QDir::currentPath(),
                                                QString("YAML files (*.yml *.yaml *.json)")).toStdString();
    }

    auto fileString = filename_;
    // cutting the path, so that only the actual filename remains
    fileString.erase(0,fileString.rfind('/') + 1);
    moleculeWidget->fileInfoText_->setText(fileString.c_str());

    YAML::Node doc = YAML::LoadAllFromFile(filename_)[1]; // load results
    auto atoms = doc["Atoms"].as<AtomsVector>();

    /* Note on the compatibility mode:
     * Clusters from result files produced without initial electron indices shuffling contain
     * artificial correlation between the core electrons, presumably originating from the Hungarian selecting the first
     * viable permutation. These artificial correlations are removed manually in the visualization.
     */
    if(doc["CompatabilityMode"] && doc["CompatabilityMode"].as<bool>())
        moleculeWidget->activateCompatabilityMode();

    /*
    auto nElectrons = doc["Clusters"][0]["Structures"][0].as<ElectronsVector>().numberOfEntities();

    auto EnStats = doc["En"].as<VectorStatistics>();
    maximaProcessingWidget->setAtomEnergies(EnStats);
    maximaProcessingWidget->setAtomsVector(atoms);
    maximaProcessingWidget->initializeTreeItems(maximaProcessingWidget->atomsTreeWidget(), int(atoms.numberOfEntities()));
    maximaProcessingWidget->initializeTreeItems(maximaProcessingWidget->electronsTreeWidget(), int(nElectrons));
    */

    moleculeWidget->setSharedAtomsVector(atoms);

    // load camera settings
    if(doc[Camera::settings.name()]) {
        Camera::settings = Settings::Camera(doc);
    }

    for (int clusterId = 0; clusterId < static_cast<int>(doc["Clusters"].size()); ++clusterId) {
        spdlog::info("{} out of {} clusters loaded...", clusterId+1, static_cast<int>(doc["Clusters"].size()));

        clusterCollection_.emplace_back(doc["Clusters"][clusterId].as<ClusterData>());
        const auto & cluster = clusterCollection_.back();

        auto item = new IntegerSortedTreeWidgetItem(
                maximaList, {QString::number(clusterId),
                 QString::number(1.0 * cluster.N_ / doc["NSamples"].as<unsigned>(), 'f', 4),
                 QString::number(cluster.valueStats_.cwiseMin()[0]/2.0, 'f', 3),
                 QString::number(cluster.valueStats_.cwiseMax()[0]/2.0, 'f', 3)});

        item->setCheckState(0, Qt::CheckState::Unchecked);

        QList<QVariant> id = {clusterId, -1};
        item->setData(0, Qt::ItemDataRole::UserRole, id);

        auto structures = doc["Clusters"][clusterId]["Structures"];

        if (structures.size() > 1) {
            for (int structureId = 0; structureId < static_cast<int>(structures.size()); ++structureId) {
                auto subItem = new IntegerSortedTreeWidgetItem(item, QStringList({
                                                                                         QString::number(structureId),
                                                                                         QString::number(1.0 *
                                                                                                         cluster.subN_[structureId] /
                                                                                                         cluster.N_,
                                                                                                         'f', 4)}));
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
}

double InPsightsWidget::sumProbabilities(){
    auto* root = maximaList->invisibleRootItem();

    double summedProbability = 0.0;

    // iterate over topLevelItems
    for (int i = 0; i < root->childCount(); ++i) {
        if(root->child(i)->checkState(0) == Qt::Checked) {
            summedProbability += root->child(i)->text(1).toDouble();
        }
        else{
            for (int j = 0; j<root->child(i)->childCount(); ++j){
                if(root->child(i)->child(j)->checkState(0) == Qt::Checked){
                    summedProbability += root->child(i)->child(j)->text(1).toDouble()*
                            root->child(i)->text(1).toDouble();
                }
            }
        }
    }
    return summedProbability;
}

void InPsightsWidget::initialView() {
    spinCorrelationsCheckBox->setTristate();
    spinCorrelationsCheckBox->setCheckState(Qt::CheckState::PartiallyChecked);
    maximaList->sortItems(0,Qt::SortOrder::AscendingOrder);
    atomsCheckBox->setCheckState(Qt::CheckState::Checked);
    bondsCheckBox->setCheckState(Qt::CheckState::Checked);
    axesCheckBox->setCheckState(Qt::CheckState::Unchecked);
    maximaList->topLevelItem(0)->setCheckState(0, Qt::CheckState::Checked);
    moleculeWidget->initialCameraSetup(
            Camera::settings.zoom.get(),
            Camera::settings.pan.get(),
            Camera::settings.tilt.get(),
            Camera::settings.roll.get()
            );
}

std::string InPsightsWidget::filenameWithoutExtension(){
    size_t lastindex = filename_.find_last_of(".");
    return filename_.substr(0, lastindex);
};

bool InPsightsWidget::plotAllActiveQ() {
    return plotAllCheckBox->checkState() == Qt::Checked;
}

void InPsightsWidget::onSedsExport(bool) {
    for (auto &sed : moleculeWidget->activeSedsMap_){
        auto clusterId = sed.first;
        int i = 0;
        for (auto &voxelCube : clusterCollection_[clusterId].voxelCubes_){
            std::string fileName = "sed-"+std::to_string(clusterId)+"-"+std::to_string(i)+".txt";
            std::string comment = filename_ + " cluster " + std::to_string(clusterId) + " SED " + std::to_string(i);
            voxelCube.exportMacmolplt(fileName, comment);
            i += 1;
        }
    }
    spdlog::info("Exported SEDs");
}

void InPsightsWidget::onDeselectAll(bool) {
    auto* root = maximaList->invisibleRootItem();
    // iterate over topLevelItems
    for (int i = 0; i < root->childCount(); ++i) {
        root->child(i)->setCheckState(0, Qt::Unchecked);
        for (int j = 0; j<root->child(i)->childCount(); ++j){
            root->child(i)->child(j)->setCheckState(0, Qt::Unchecked);
        }
    }
}
