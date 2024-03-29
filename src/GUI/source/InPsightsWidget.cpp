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
        atomsCheckBox(new QCheckBox("Nuclei", this)),
        bondsCheckBox(new QCheckBox("Bonds", this)),
        axesCheckBox(new QCheckBox("Axes", this)),
        sampleAverageCheckBox(new QCheckBox("Sample average", this)),
        spinCorrelationsCheckBox(new QCheckBox("Spin correlations", this)),
        sedsCheckBox(new QCheckBox("SEDs", this)),
        maximaHullsCheckBox(new QCheckBox("Maxima hulls", this)),
        plotAllCheckBox(new QCheckBox("All of cluster", this)),
        coloredCheckBox(new QCheckBox("Multicolored", this)),
        moveElectronsCheckBox(new QCheckBox("Move electrons", this)),
        electronsNumberCheckBox(new QCheckBox("Display indices", this)),
        spinCorrelationBox(new QDoubleSpinBox(this)),
        sedPercentageBox(new QDoubleSpinBox(this)),
        scaleVectorBox(new QDoubleSpinBox(this)),
        bondBox(new QDoubleSpinBox(this)),
        atom1Box(new QSpinBox(this)),
        atom2Box(new QSpinBox(this)),
        eigenvectorSpinBox(new QSpinBox(this)),
        electron1Box(new QSpinBox(this)),
        electron2Box(new QSpinBox(this)),
        maximaList(new QTreeWidget(this)),
        probabilitySum(new QLabel(this)),
        eigenvalueLabel(new QLabel(this)),
        deselectAllButton(new QPushButton("Deselect all", this)),
        lastMovedElectronClusterVector({0, 0, -1}),
        globalMinPhi(0.0f)
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
    auto vboxSettings = new QVBoxLayout();
    auto vboxParticleHighlighting = new QVBoxLayout();
    auto vboxEigenvectors = new QVBoxLayout();
    auto selectorBox = new QGroupBox("Particle highlighting:");
    auto settingsBox = new QGroupBox("Settings:");
    auto eigenvectorsBox = new QGroupBox("Eigenvectors:");

    setLayout(hbox);

    resize(1424, 801);
    hbox->addWidget(moleculeWidget, 1);
    hbox->addLayout(vboxOuter, 0);

    maximaList->setMinimumWidth(352);

    // put into MaximaTreeWidget class
    auto headerLabels = QList<QString>({"ID", "Weight", "min ΔΦ", "max ΔΦ"});
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
    selectorBox->setLayout(vboxParticleHighlighting);

    auto selectorGrid = new QGridLayout();
    vboxParticleHighlighting->addLayout(selectorGrid,1);
    selectorGrid->addWidget(new QLabel("Nuclei"),0,0);
    selectorGrid->addWidget(atom1Box,0,1);
    selectorGrid->addWidget(atom2Box,0,2);
    selectorGrid->addWidget(new QLabel("Electrons"),1,0);
    selectorGrid->addWidget(electron1Box,1,1);
    selectorGrid->addWidget(electron2Box,1,2);

    if (not clusterCollection_[0].eigenvalues_.empty()) {
        vboxOuter->addWidget(eigenvectorsBox);
        eigenvectorsBox->setLayout(vboxEigenvectors);

        auto eigenvectorsGrid = new QGridLayout();
        vboxEigenvectors->addLayout(eigenvectorsGrid, 1);

        // first row
        eigenvectorsGrid->addWidget(eigenvectorSpinBox,0,0);
        eigenvectorsGrid->addWidget(eigenvalueLabel,0,1);
        resetEigenvalueLabel();

        //second row
        eigenvectorsGrid->addWidget(moveElectronsCheckBox,1,1);
        eigenvectorsGrid->addWidget(scaleVectorBox,1,0);
    }
    else{
        eigenvectorSpinBox->deleteLater();
        eigenvalueLabel->deleteLater();
        moveElectronsCheckBox->deleteLater();
        scaleVectorBox->deleteLater();
    }
    
    //vboxOuter->addWidget(maximaProcessingWidget,1);
    vboxOuter->addWidget(settingsBox);
    settingsBox->setLayout(vboxSettings);

    maximaList->setSortingEnabled(true);

    auto checkboxGrid = new QGridLayout();
    vboxSettings->addLayout(checkboxGrid,1);

    // first row
    checkboxGrid->addWidget(atomsCheckBox,0,0);
    checkboxGrid->addWidget(axesCheckBox,0,1);
    axesCheckBox->setToolTip("x, y, and z axes in red, green, and blue, respectively.");

    //second row
    checkboxGrid->addWidget(bondsCheckBox,1,0);
    checkboxGrid->addWidget(bondBox,1,1);

    //third row
    checkboxGrid->addWidget(sampleAverageCheckBox,2,0);
    checkboxGrid->addWidget(coloredCheckBox,2,1);


    //fourth row
    checkboxGrid->addWidget(maximaHullsCheckBox,3,0);
    checkboxGrid->addWidget(plotAllCheckBox,3,1);

    //fifth row
    checkboxGrid->addWidget(sedsCheckBox,4,0);
    checkboxGrid->addWidget(sedPercentageBox,4,1);

    //sixth row
    checkboxGrid->addWidget(spinCorrelationsCheckBox,5,0);
    checkboxGrid->addWidget(spinCorrelationBox,5,1);

    //seventh row
    checkboxGrid->addWidget(electronsNumberCheckBox,6,0);

    setupSpinBoxes();
    setupLabels();
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

    connect(plotAllCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onPlotAllChecked);

    connect(spinCorrelationBox, qOverload<double>(&QDoubleSpinBox::valueChanged),
            this, &InPsightsWidget::onSpinCorrelationsBoxChanged);

    connect(electronsNumberCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onIndicesChecked);

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

    connect(sedPercentageBox,qOverload<double>(&QDoubleSpinBox::valueChanged),
            this, &InPsightsWidget::onSedBoxChanged);

    connect(bondBox,qOverload<double>(&QDoubleSpinBox::valueChanged),
            this, &InPsightsWidget::onBondBoxChanged);

    connect(coloredCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::updateSelectedStructures);

    connect(maximaHullsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::updateSelectedStructures);

    connect(sampleAverageCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onSampleAverageCheckBoxChanged);

    connect(eigenvectorSpinBox, qOverload<int>(&QSpinBox::valueChanged),
            this, &InPsightsWidget::onEigenvectorSpinBoxChanged);

    connect(scaleVectorBox, qOverload<double>(&QDoubleSpinBox::valueChanged),
            this, &InPsightsWidget::onScaleVectorBoxChanged);

    connect(moveElectronsCheckBox, &QCheckBox::stateChanged,
            this, &InPsightsWidget::onMoveElectronsCheckBoxChecked);

    connect(deselectAllButton, &QPushButton::clicked,
            this, &InPsightsWidget::onDeselectAll);
}

void InPsightsWidget::setupLabels() {
    eigenvalueLabel->setAlignment(Qt::AlignCenter);
}

void InPsightsWidget::setupSpinBoxes() {

    scaleVectorBox->setRange(-10,10);
    scaleVectorBox->setWrapping(false);
    scaleVectorBox->setSingleStep(0.1);
    scaleVectorBox->setValue(1);
    scaleVectorBox->setAccelerated(true);

    int numberEigenvalues;
    if(clusterCollection_[0].eigenvalues_.empty()) {
        numberEigenvalues = 1;
    }
    else {
        numberEigenvalues = clusterCollection_[0].eigenvalues_[0].size();
    }
    eigenvectorSpinBox->setRange(-1,numberEigenvalues-1);
    eigenvectorSpinBox->setWrapping(true);
    eigenvectorSpinBox->setSingleStep(1);
    eigenvectorSpinBox->setValue(-1);
    eigenvectorSpinBox->setSpecialValueText(tr("none"));
    eigenvectorSpinBox->setAccelerated(true);

    spinCorrelationBox->setRange(0.0,1.0);
    spinCorrelationBox->setSingleStep(0.01);
    spinCorrelationBox->setValue(1.0);
    spinCorrelationBox->setButtonSymbols(QAbstractSpinBox::PlusMinus);
    spinCorrelationBox->setAccelerated(true);

    sedPercentageBox->setRange(0.0,1.0);
    sedPercentageBox->setSingleStep(0.01);
    sedPercentageBox->setValue(0.2);
    sedPercentageBox->setButtonSymbols(QAbstractSpinBox::PlusMinus);
    sedPercentageBox->setAccelerated(true);

    bondBox->setRange(0.0,2.0);
    bondBox->setSingleStep(0.01);
    bondBox->setValue(0.5);
    bondBox->setButtonSymbols(QAbstractSpinBox::PlusMinus);
    bondBox->setAccelerated(true);
    bondBox->setToolTip("Threshold as multiple of the summed vdW radii.");

    atom1Box->setRange(-1,moleculeWidget->getAtomsNumber()-1);
    atom1Box->setWrapping(true);
    atom1Box->setSingleStep(1);
    atom1Box->setValue(-1);
    atom1Box->setSpecialValueText(tr("none"));
    atom1Box->setAccelerated(true);

    atom2Box->setRange(-1,moleculeWidget->getAtomsNumber()-1);
    atom2Box->setWrapping(true);
    atom2Box->setSingleStep(1);
    atom2Box->setValue(-1);
    atom2Box->setSpecialValueText(tr("none"));
    atom2Box->setAccelerated(true);

    electron1Box->setRange(-1,clusterCollection_[0].representativeStructure().numberOfEntities()-1);
    electron1Box->setWrapping(true);
    electron1Box->setSingleStep(1);
    electron1Box->setValue(-1);
    electron1Box->setSpecialValueText(tr("none"));
    electron1Box->setAccelerated(true);

    electron2Box->setRange(-1,clusterCollection_[0].representativeStructure().numberOfEntities()-1);
    electron2Box->setWrapping(true);
    electron2Box->setSingleStep(1);
    electron2Box->setValue(-1);
    electron2Box->setSpecialValueText(tr("none"));
    electron2Box->setAccelerated(true);
}

void InPsightsWidget::selectedStructure(QTreeWidgetItem *item, int column) {
    //auto maximaTreeWidgetItem = dynamic_cast<IntegerSortedTreeWidgetItem*>(item);

    if (column != 0)
        spdlog::critical("Column 0 expected but got {} ", column);
    std::vector<int> TickedStructuresCountVector = getTickedStructuresCountVector();
    int count = TickedStructuresCountVector[0];
    bool checkEigenval = checkEigenvalues();
    auto id = item->data(0, Qt::ItemDataRole::UserRole).toList();
    auto clusterId = id[0].toInt();
    auto secondId = id[1].toInt();
    auto structureId = secondId;
    if (structureId == -1){  // if secondId == -1, the top level checkbox is marked
        structureId = 0;  // taking structure 0 to visualize the top level cluster
    }

    auto createQ = item->checkState(0) == Qt::CheckState::Checked;

    if (createQ) {
        if (not clusterCollection_[0].eigenvalues_.empty()) {
            if (eigenvectorSpinBox->value() != -1 and count > 1) {
                spdlog::warn("Chose only one structure for eigenvalues");
            }
        }
        if (sampleAverageCheckBox->checkState() == Qt::CheckState::Checked){
            if (not clusterCollection_[0].eigenvalues_.empty()) {
                if (eigenvectorSpinBox->value() != -1) {
                    spdlog::warn("No eigenvectors were calculated for the sample average");
                    eigenvectorSpinBox->setValue(-1);
                }
                if (moveElectronsCheckBox->isChecked()) {
                    spdlog::warn("Cannot move electrons for the sample average");
                    moveElectronsCheckBox->setCheckState(Qt::Unchecked);
                }
            }
            auto sampleAverage = clusterCollection_[clusterId].sampleAverage_;
            if (sampleAverage.numberOfEntities() > 0) {
                moleculeWidget->addElectronsVector(sampleAverage, clusterId, secondId,
                                                   coloredCheckBox->checkState() == Qt::Checked);
            }
            else {
                moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId],
                                                   clusterId, secondId, coloredCheckBox->checkState() == Qt::Checked);
                spdlog::warn("Sample averaged vectors were not calculated. Plotting the first maximum instead...");
            }
        }
        else {
            moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId],
                                               clusterId, secondId, coloredCheckBox->checkState() == Qt::Checked);
        }
        if (not clusterCollection_[0].eigenvalues_.empty()) {
            if (checkEigenvalues() and eigenvectorSpinBox->value() != -1) {
                eigenvalueLabel->setText(QString::number(clusterCollection_[clusterId].eigenvalues_[structureId][eigenvectorSpinBox->value()] / 2, 'f', 4));
                moleculeWidget->drawEigenvectors(true, clusterId, structureId, eigenvectorSpinBox->value(), scaleVectorBox->value());
            }
            else {
                if (moveElectronsCheckBox->isChecked()) {
                    moveElectronsCheckBox->setCheckState(Qt::Unchecked);
                    std::cout << lastMovedElectronClusterVector[0] << " " << lastMovedElectronClusterVector[1] << " " << lastMovedElectronClusterVector[2] << std::endl;
                    moleculeWidget->removeElectronsVector(lastMovedElectronClusterVector[0], lastMovedElectronClusterVector[2]);
                    moleculeWidget->addElectronsVector(clusterCollection_[lastMovedElectronClusterVector[0]]
                                                        .exemplaricStructures_[lastMovedElectronClusterVector[1]],
                                                       lastMovedElectronClusterVector[0], lastMovedElectronClusterVector[2],
                                                       coloredCheckBox->checkState() == Qt::Checked);
                }
                moleculeWidget->removeEigenvectors();
                eigenvalueLabel->setText(QString(" "));
            }
        }
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
    }
    else {
        if (not clusterCollection_[0].eigenvalues_.empty()) {
            if (eigenvectorSpinBox->value() != -1 and count > 1) {
                spdlog::warn("Chose only one structure for eigenvalues");
            }
            if (not checkEigenval) {
                moveElectronsCheckBox->setCheckState(Qt::Unchecked);
            }
            if (checkEigenval and eigenvectorSpinBox->value() != -1) {
                auto* root = maximaList->invisibleRootItem();
                unsigned cluId = -1;
                unsigned strId = 0;
                // iterate over topLevelItems
                for (int i = 0; i < root->childCount(); ++i) {
                    if(root->child(i)->checkState(0) == Qt::Checked) {
                        cluId = i;
                    }
                    // iterate over childs of topLevelItem i
                    for (int j = 0; j<root->child(i)->childCount(); ++j) {
                        if(root->child(i)->child(j)->checkState(0) == Qt::Checked) {
                            cluId = i;
                            strId = j;
                        }
                    }
                }
                eigenvalueLabel->setText(QString::number(clusterCollection_[cluId].eigenvalues_[strId][eigenvectorSpinBox->value()] / 2, 'f', 4));
                moleculeWidget->drawEigenvectors(true, cluId, strId, eigenvectorSpinBox->value(), scaleVectorBox->value());
            }
            else {
                moleculeWidget->removeEigenvectors();
                eigenvalueLabel->setText(QString(" "));
            }
        }

        moleculeWidget->removeElectronsVector(clusterId, secondId);

        if (moleculeWidget->activeSedsMap_.find(clusterId) != moleculeWidget->activeSedsMap_.end())
            moleculeWidget->removeSeds(clusterId);

        if (moleculeWidget->activeMaximaHullsMap_.find(clusterId) != moleculeWidget->activeMaximaHullsMap_.end())
            moleculeWidget->removeMaximaHulls(clusterId);
    }
    redrawSpinDecorations();
    probabilitySum->setText(QString("Σ Weight = ") + QString::number(sumProbabilities(), 'f', 4));
};

void InPsightsWidget::resetEigenvalueLabel() {
    eigenvalueLabel->setText(QString(" "));
}

void InPsightsWidget::addMovedElectronsVector(int clusterId, int structureId, int secondId) {
    unsigned electronsNumber = clusterCollection_[clusterId].exemplaricStructures_[structureId].numberOfEntities();
    auto startElectronsVector = clusterCollection_[clusterId].exemplaricStructures_[structureId];
    auto eigenvector = clusterCollection_[clusterId].eigenvectors_[structureId];
    PositionsVector movedElectronPositions;
    for (unsigned i = 0; i < electronsNumber; ++i) {
        movedElectronPositions.append(startElectronsVector.positionsVector()[i] +
                                      eigenvector[eigenvectorSpinBox->value() * electronsNumber + i] *
                                      scaleVectorBox->value());
    }
    auto movedElectronsVector = ElectronsVector(movedElectronPositions,
                                                startElectronsVector.typesVector());
    moleculeWidget->addElectronsVector(movedElectronsVector, clusterId, secondId);
    lastMovedElectronClusterVector[0] = clusterId;
    lastMovedElectronClusterVector[1] = structureId;
    lastMovedElectronClusterVector[2] = secondId;
}

bool InPsightsWidget::checkEigenvalues() {
    if(clusterCollection_[0].eigenvalues_.empty()) {
        return false;
    }
    else {
        std::vector<int> TickedStructuresCountVector = getTickedStructuresCountVector();
        int count = TickedStructuresCountVector[0];
        if (count > 1) {
            return false;
        }
        if (count == 0) {
            return false;
        }
    }
    return true;
}

std::vector<int> InPsightsWidget::getTickedStructuresCountVector() {
    auto* root = maximaList->invisibleRootItem();
    std::vector<int> TickedStructuresCountVector;
    unsigned count = 0;
    unsigned clusterId = -1;
    unsigned structureId = -1;
    // iterate over topLevelItems
    for (int i = 0; i < root->childCount(); ++i) {
        if(root->child(i)->checkState(0) == Qt::Checked) {
            count += 1;
            clusterId = i;
        }
        // iterate over childs of topLevelItem i
        for (int j = 0; j<root->child(i)->childCount(); ++j) {
            if(root->child(i)->child(j)->checkState(0) == Qt::Checked) {
                count += 1;
                clusterId = i;
                structureId = j;
            }
        }
    }
    TickedStructuresCountVector.emplace_back(count);
    TickedStructuresCountVector.emplace_back(clusterId);
    TickedStructuresCountVector.emplace_back(structureId);
    return TickedStructuresCountVector;
}

void InPsightsWidget::onSampleAverageCheckBoxChanged(int stateId) {
    if (not clusterCollection_[0].eigenvalues_.empty()) {
        if (moveElectronsCheckBox->isChecked()) {
            sampleAverageCheckBox->setCheckState(Qt::Unchecked);
            spdlog::warn("Cannot move electrons for the sample average");
        }
        if (eigenvectorSpinBox->value() != -1) {
            sampleAverageCheckBox->setCheckState(Qt::Unchecked);
            spdlog::warn("No eigenvectors were calculated for the sample average");
        }
        if (not moveElectronsCheckBox->isChecked() and eigenvectorSpinBox->value() == -1) {
            updateSelectedStructures(42);
        }
    }
    else {
        updateSelectedStructures(42);
    }

}

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

void InPsightsWidget::onIndicesChecked(int stateId) {
    moleculeWidget->drawIndices(Qt::CheckState(stateId) == Qt::CheckState::Checked);
}

void InPsightsWidget::onBondsChecked(int stateId) {
    if (atomsCheckBox->isChecked()) {
        moleculeWidget->drawBonds(Qt::CheckState(stateId) == Qt::CheckState::Checked, bondBox->value());
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
    if (clusterCollection_[0].SeeStats_.getTotalWeight() == 0) {
        if ((Qt::CheckState(stateId) == Qt::CheckState::Checked) or (Qt::CheckState(stateId) == Qt::CheckState::PartiallyChecked)) {
            spdlog::warn("Spin correlations were not calculated.");
            spinCorrelationsCheckBox->setCheckState(Qt::CheckState::Unchecked);
        }
    }
    else {
        moleculeWidget->deleteSpinCorrelations();
        if(Qt::CheckState(stateId) == Qt::Checked)
            moleculeWidget->drawSpinCorrelations(clusterCollection_, spinCorrelationBox->value(),true);
        else if(Qt::CheckState(stateId) == Qt::PartiallyChecked)
            moleculeWidget->drawSpinCorrelations(clusterCollection_, spinCorrelationBox->value(),false);
    }
}

void InPsightsWidget::onSpinCorrelationsBoxChanged(double value){
    if (clusterCollection_[0].SeeStats_.getTotalWeight() > 0)
        redrawSpinDecorations();
}

void InPsightsWidget::onBondBoxChanged(double value) {
    if (bondsCheckBox->isChecked()) {
        bondsCheckBox->setCheckState(Qt::Unchecked);
        bondsCheckBox->setCheckState(Qt::Checked);
    }
}

void InPsightsWidget::onAtom1BoxChanged(int value) {
    if (atomsCheckBox->isChecked()) {
        moleculeWidget->onAtomsHighlighted(value);
    }
}

void InPsightsWidget::onAtom2BoxChanged(int value) {
    if (atomsCheckBox->isChecked()) {
        moleculeWidget->onAtomsChecked(value);
    }
}

void InPsightsWidget::onElectron1BoxChanged(int value) {
    moleculeWidget->onElectronsHighlighted(value);
}

void InPsightsWidget::onElectron2BoxChanged(int value) {
    moleculeWidget->onElectronsChecked(value);
}

void InPsightsWidget::onEigenvectorSpinBoxChanged(int value) {
    if (sampleAverageCheckBox->isChecked() and value != -1) {
        eigenvectorSpinBox->setValue(-1);
        value = -1;
        spdlog::warn("No eigenvectors were calculated for the sample average");
    }
    std::vector<int> tickedStructuresCountVector = getTickedStructuresCountVector();
    int count = tickedStructuresCountVector[0];
    int id = tickedStructuresCountVector[1];
    int secondId = tickedStructuresCountVector[2];
    int structureId = secondId;
    if (structureId == -1) {
        structureId = 0;
    }
    if (value == -1) {
        if (not sampleAverageCheckBox->isChecked() and moveElectronsCheckBox->isChecked()) {
            moveElectronsCheckBox->setCheckState(Qt::Unchecked);
            spdlog::warn("Chose an eigenvector to move the electrons");
            moleculeWidget->removeElectronsVector(id, secondId);
            moleculeWidget->addElectronsVector(clusterCollection_[id].exemplaricStructures_[structureId],
                                               id, secondId,
                                               coloredCheckBox->checkState() == Qt::Checked);
        }
    }
    if (value != -1) {
        if (count > 1) {
            spdlog::warn("Chose only one structure for eigenvalues");
        }
        else if (count == 1){
            eigenvalueLabel->setText(QString::number(clusterCollection_[id].eigenvalues_[structureId][value] / 2, 'f', 4));
            moleculeWidget->drawEigenvectors(true, id, structureId, value, scaleVectorBox->value());
            if (moveElectronsCheckBox->isChecked()) {
                moveElectronsCheckBox->setCheckState(Qt::Unchecked);
            }
        }
        if (count == 0) {
            eigenvalueLabel->setText(QString(" "));
            spdlog::warn("Chose a structure for eigenvalues");
        }
    }
    if (value == -1 ) {
        eigenvalueLabel->setText(QString(" "));
        moleculeWidget->removeEigenvectors();
    }
}

void InPsightsWidget::onScaleVectorBoxChanged(double value) {
    std::vector<int> tickedStructuresCountVector = getTickedStructuresCountVector();
    int count = tickedStructuresCountVector[0];
    if (count == 1 and eigenvectorSpinBox->value() != -1) {
        int clusterId = tickedStructuresCountVector[1];
        int structureId = tickedStructuresCountVector[2];
        int secondId = structureId;
        if (structureId == -1) {
            structureId = 0;
        }
        if (moveElectronsCheckBox->checkState() == Qt::Checked) {
            moleculeWidget->removeElectronsVector(clusterId, secondId);
            addMovedElectronsVector(clusterId, structureId, secondId);
            moleculeWidget->drawEigenvectors(true, clusterId, structureId, eigenvectorSpinBox->value(), value);
        }
        else {
            moleculeWidget->drawEigenvectors(true, clusterId, structureId, eigenvectorSpinBox->value(), value);
        }
    }
    redrawSpinDecorations();
}

void InPsightsWidget::onMoveElectronsCheckBoxChecked(int stateId){
    if (sampleAverageCheckBox->isChecked() and moveElectronsCheckBox->isChecked()) {
        moveElectronsCheckBox->setCheckState(Qt::Unchecked);
        spdlog::warn("Cannot move electrons for the sample average");
    }
    std::vector<int> tickedStructuresCountVector = getTickedStructuresCountVector();
    int count = tickedStructuresCountVector[0];
    if (moveElectronsCheckBox->isChecked()) {
        if (count == 0){
            spdlog::warn("Chose a structure to move the electrons");
            moveElectronsCheckBox->setCheckState(Qt::Unchecked);
        }
        if (count > 1) {
            spdlog::warn("Chose only one structure to move the electrons");
            moveElectronsCheckBox->setCheckState(Qt::Unchecked);
        }
        if (count == 1 and eigenvectorSpinBox->value() == -1) {
            spdlog::warn("Chose an eigenvector to move the electrons");
            moveElectronsCheckBox->setCheckState(Qt::Unchecked);
        }
    }

    if (count == 1 and eigenvectorSpinBox->value() != -1) {
        int clusterId = tickedStructuresCountVector[1];
        int structureId = tickedStructuresCountVector[2];
        int secondId = structureId;
        if (structureId == -1) {
            structureId = 0;
        }
        if (moveElectronsCheckBox->checkState() == Qt::Unchecked) {
            moleculeWidget->removeElectronsVector(clusterId, secondId);
            moleculeWidget->addElectronsVector(clusterCollection_[clusterId].exemplaricStructures_[structureId],
                                               clusterId, secondId,
                                               coloredCheckBox->checkState() == Qt::Checked);
        }
        else {
            moleculeWidget->removeElectronsVector(clusterId, secondId);
            addMovedElectronsVector(clusterId, structureId, secondId);
        }
    }
    redrawSpinDecorations();
}

void InPsightsWidget::showSplashScreen() {
    auto splashScreen = new QSplashScreen();
    auto pixmap = QPixmap(":inPsights.png").scaledToWidth(400, Qt::TransformationMode::SmoothTransformation);

    splashScreen->setPixmap(pixmap);
    splashScreen->show();

    std::string message = inPsights::version() + "\n"\
                          "LuFG Theoretical Chemistry, RWTH Aachen University";

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
     * artificial correlation between the core (or ionic) electrons, presumably originating from the
     * Hungarian selecting the first viable permutation. These artificial correlations are removed
     * manually in the visualization.
     */
    if(!doc["CompatibilityMode"] or doc["CompatibilityMode"].as<bool>()) {
        moleculeWidget->activateCompatibilityMode();
    }

    moleculeWidget->setSharedAtomsVector(atoms);

    // load camera settings
    if(doc[Camera::settings.name()]) {
        Camera::settings = Settings::Camera(doc);
    }

    if (doc["GlobalMinPhi"]) {
        globalMinPhi = doc["GlobalMinPhi"].as<float>();
    }
    else {
        globalMinPhi = doc["Clusters"][0].as<ClusterData>().valueStats_.cwiseMin()[0]/2.0;
    }
    for (int clusterId = 0; clusterId < static_cast<int>(doc["Clusters"].size()); ++clusterId) {
        spdlog::info("{} out of {} clusters loaded...", clusterId+1, static_cast<int>(doc["Clusters"].size()));

        clusterCollection_.emplace_back(doc["Clusters"][clusterId].as<ClusterData>());
        const auto & cluster = clusterCollection_.back();
        auto phiStrings = getPhiStrings(cluster.valueStats_);

        auto item = new IntegerSortedTreeWidgetItem(
                maximaList, {QString::number(clusterId),
                 QString::number(1.0 * cluster.N_ / doc["NSamples"].as<unsigned>(), 'f', 4),
                 phiStrings.first,
                 phiStrings.second});

        item->setCheckState(0, Qt::CheckState::Unchecked);

        QList<QVariant> id = {clusterId, -1};
        item->setData(0, Qt::ItemDataRole::UserRole, id);

        auto structures = doc["Clusters"][clusterId]["Structures"];

        if (structures.size() > 1) {
            for (int structureId = 0; structureId < static_cast<int>(structures.size()); ++structureId) {
                if (not cluster.subValueStats_.empty()){
                    phiStrings = getPhiStrings(cluster.subValueStats_[structureId]);
                }else{
                    phiStrings = std::pair<QString, QString>(QString(' '),QString(' '));
                }

                auto subItem = new IntegerSortedTreeWidgetItem(item,
                                                               QStringList({
                                                                  QString::number(structureId),
                                                                  QString::number(1.0 *
                                                                                  cluster.subN_[structureId] /
                                                                                  cluster.N_,
                                                                                  'f', 4),
                                                                  phiStrings.first,
                                                                  phiStrings.second}));
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

std::pair<QString, QString> InPsightsWidget::getPhiStrings(const SingleValueStatistics& valueStats){
    float minPhi = valueStats.cwiseMin()[0]/2.0;
    auto minPhiString = QString::number(minPhi-globalMinPhi, 'f', 4);
    if (minPhi-globalMinPhi >= 0)
        minPhiString = QString(' ') + minPhiString;

    float maxPhi = valueStats.cwiseMax()[0]/2.0;
    auto maxPhiString = QString::number(maxPhi-globalMinPhi, 'f', 4);
    if (maxPhi-globalMinPhi >= 0)
        maxPhiString = QString(' ') + maxPhiString;
    return std::pair<QString, QString>(minPhiString.left(7), maxPhiString.left(7));
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

void InPsightsWidget::onSedBoxChanged(double) {
    if (sedsCheckBox->isChecked())
        updateSelectedStructures(Qt::CheckState::Checked);
}
