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
#include <Version.h>
#include <vector>
#include <ParticlesVector.h>
#include "CameraSettings.h"
#include <spdlog/spdlog.h>


InPsightsWidget::InPsightsWidget(QWidget *parent, const std::string& filename)
        :
        QWidget(parent),
        filename_(filename),
        moleculeWidget(new MoleculeWidget(this)),
        maximaProcessingWidget(new MaximaProcessingWidget(this)),
        atomsCheckBox(new QCheckBox("Atoms", this)),
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
    auto gbox = new QGroupBox("Settings:");

    setLayout(hbox);

    resize(1280, 800);
    hbox->addWidget(moleculeWidget, 2);
    hbox->addLayout(vboxOuter, 1);

    // put into MaximaTreeWidget class
    auto headerLabels = QList<QString>({"ID", "P", "min(Φ)", "max(Φ)"});
    maximaList->setColumnCount(headerLabels.size());
    maximaList->setHeaderLabels(headerLabels);
    maximaList->header()->setStretchLastSection(false);

    maximaList->headerItem()->setToolTip(0,QString("IDs sorted by Φ value"));
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

    vboxOuter->addWidget(maximaProcessingWidget,1);
    vboxOuter->addWidget(gbox);
    gbox->setLayout(vboxInner);

    maximaList->setSortingEnabled(true);

    auto checkboxGrid = new QGridLayout();
    vboxInner->addLayout(checkboxGrid,1);

    // first column
    checkboxGrid->addWidget(atomsCheckBox,0,0);
    checkboxGrid->addWidget(bondsCheckBox,1,0);
    checkboxGrid->addWidget(axesCheckBox,2,0);
    checkboxGrid->addWidget(sedsExportButton,3,0);

    // second column
    checkboxGrid->addWidget(sampleAverageCheckBox,0,1);
    checkboxGrid->addWidget(maximaHullsCheckBox,1,1);
    checkboxGrid->addWidget(spinCorrelationsCheckBox,2,1);
    checkboxGrid->addWidget(sedsCheckBox,3,1);

    // third column
    checkboxGrid->addWidget(plotAllCheckBox,0,2);
    checkboxGrid->addWidget(coloredCheckBox,1,2);
    checkboxGrid->addWidget(spinCorrelationBox,2,2);
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

    connect(spinCorrelationsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSpinCorrelationsChecked);

    connect(spinCorrelationsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSpinCorrelationsChecked);

    connect(plotAllCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onPlotAllChecked);


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
}

void InPsightsWidget::selectedStructure(QTreeWidgetItem *item, int column) {
    //auto maximaTreeWidgetItem = dynamic_cast<IntegerSortedTreeWidgetItem*>(item);

    if (column != 0)
        spdlog::critical("Column 0 expected but got {} ", column);

    auto id = item->data(0, Qt::ItemDataRole::UserRole).toList();
    auto clusterId = id[0].toInt();
    auto structureId = id[1].toInt();

    auto createQ = item->checkState(0) == Qt::CheckState::Checked;

    if (createQ) {
        auto sampleAverage = clusterCollection_[clusterId].sampleAverage_;

        if(sampleAverageCheckBox->checkState() == Qt::CheckState::Checked && sampleAverage.numberOfEntities() > 0) {
            moleculeWidget->addElectronsVector(sampleAverage, clusterId, structureId, coloredCheckBox->checkState() == Qt::Checked);
        } else if (sampleAverageCheckBox->checkState() == Qt::CheckState::Checked
        && sampleAverage.numberOfEntities() == 0) {
            moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId],
                                               clusterId, structureId, coloredCheckBox->checkState() == Qt::Checked);
            spdlog::warn("Sample averaged vectors were not calculated. Plotting the first maximum instead...");
        } else {
            moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId],
                                               clusterId, structureId, coloredCheckBox->checkState() == Qt::Checked);
        }

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
    probabilitySum->setText(QString("Σ=") + QString::number(sumProbabilities(), 'f', 4));
};

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
    if(spinCorrelationsCheckBox->checkState() == Qt::Checked){
        onSpinCorrelationsChecked(Qt::Unchecked);
        onSpinCorrelationsChecked(Qt::Checked);
    }
    else if(spinCorrelationsCheckBox->checkState() == Qt::PartiallyChecked){
        onSpinCorrelationsChecked(Qt::Unchecked);
        onSpinCorrelationsChecked(Qt::PartiallyChecked);
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

void InPsightsWidget::onSpinCorrelationsChecked(int stateId) {
    if(Qt::CheckState(stateId) == Qt::Checked)
        moleculeWidget->drawSpinCorrelations(true,clusterCollection_, spinCorrelationBox->value(),true);
    else if(Qt::CheckState(stateId) == Qt::PartiallyChecked)
        moleculeWidget->drawSpinCorrelations(true,clusterCollection_, spinCorrelationBox->value(),false);
    else
        moleculeWidget->drawSpinCorrelations(false,clusterCollection_, spinCorrelationBox->value(),false);
}

void InPsightsWidget::onSpinCorrelationsBoxChanged(double value) {
    if (spinCorrelationsCheckBox->checkState() == Qt::CheckState::Checked) {
        onSpinCorrelationsChecked(Qt::CheckState::Unchecked); //TODO ugly, create update() function in SpinCorrelation3D and make it accessible
        onSpinCorrelationsChecked(Qt::CheckState::Checked);
    } else if (spinCorrelationsCheckBox->checkState() == Qt::CheckState::PartiallyChecked) {
            onSpinCorrelationsChecked(Qt::CheckState::Unchecked);
            onSpinCorrelationsChecked(Qt::CheckState::PartiallyChecked);
    }else {
         onSpinCorrelationsChecked(Qt::CheckState::Unchecked);
    }
}

void InPsightsWidget::showSplashScreen() {
    auto splashScreen = new QSplashScreen();
    auto pixmap = QPixmap(":inPsights.png").scaledToWidth(400, Qt::TransformationMode::SmoothTransformation);

    splashScreen->setPixmap(pixmap);
    splashScreen->show();

    std::string message = inPsights::version() + " (alpha)\n "\
                          "Copyright © 2016-2020  Michael A. Heuer.";
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

    moleculeWidget->fileInfoText_->setText(filename_.c_str());

    YAML::Node doc = YAML::LoadAllFromFile(filename_)[1]; // load results
    auto atoms = doc["Atoms"].as<AtomsVector>();

    /* Note on the compatibility mode:
     * Clusters from result files produced without initial electron indices shuffling contain
     * artificial correlation between the core electrons, presumably originating from the Hungarian selecting the first
     * viable permutation. These artificial correlations are removed manually in the visualization.
     */
    if(doc["CompatabilityMode"] && doc["CompatabilityMode"].as<bool>())
        moleculeWidget->activateCompatabilityMode();


    auto nElectrons = doc["Clusters"][0]["Structures"][0].as<ElectronsVector>().numberOfEntities();

    auto EnStats = doc["En"].as<VectorStatistics>();
    maximaProcessingWidget->setAtomEnergies(EnStats);
    maximaProcessingWidget->setAtomsVector(atoms);
    maximaProcessingWidget->initializeTreeItems(maximaProcessingWidget->atomsTreeWidget(), int(atoms.numberOfEntities()));
    maximaProcessingWidget->initializeTreeItems(maximaProcessingWidget->electronsTreeWidget(), int(nElectrons));


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

        QList<QVariant> id = {clusterId, 0};
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
