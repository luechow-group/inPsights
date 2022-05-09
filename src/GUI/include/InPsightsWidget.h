// Copyright (C) 2018-2019 Michael Heuer.
// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_INPSIGHTSWIDGET_H
#define INPSIGHTS_INPSIGHTSWIDGET_H

#include <MoleculeWidget.h>

#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QTreeWidget>
#include <ClusterData.h>

class InPsightsWidget : public QWidget {
Q_OBJECT
public:
    explicit InPsightsWidget(QWidget *parent = nullptr, const std::string &filename = "");
    bool plotAllActiveQ();
    std::string filenameWithoutExtension();
    std::vector<ClusterData> clusterCollection_;

public slots:
    void selectedStructure(QTreeWidgetItem *item, int column);
    void onAtomsChecked(int stateId);
    void onBondsChecked(int stateId);
    void onAxesChecked(int stateId);
    void onSpinCorrelationsChecked(int stateId= 0);
    void onPlotAllChecked(int stateId= 0);
    void onSpinCorrelationsBoxChanged(double value);
    void onElectron1BoxChanged(int value);
    void onElectron2BoxChanged(int value);
    void onAtom1BoxChanged(int value);
    void onAtom2BoxChanged(int value);
    void onDeselectAll(bool);
    void onSedChecked(int stateId);
    void onSedBoxChanged(double value);
    void onEigenvectorSpinBoxChanged(int value);
    void onScaleVectorBoxChanged(double value);
    void onBondBoxChanged(double value);
    void onMoveElectronsCheckBoxChecked(int stateId);
    void updateSelectedStructures(int);
    void addMovedElectronsVector(int clusterId, int structureId, int secondId);
    bool checkEigenvalues();
    void onSampleAverageCheckBoxChanged(int stateId);
    std::vector<int> getTickedStructuresCountVector();

private:
    std::string filename_;
    MoleculeWidget *moleculeWidget;
    QCheckBox *atomsCheckBox, *bondsCheckBox, *axesCheckBox, *sampleAverageCheckBox, *spinCorrelationsCheckBox,
    *sedsCheckBox,*maximaHullsCheckBox, *plotAllCheckBox, *coloredCheckBox, *moveElectronsCheckBox;
    QDoubleSpinBox *spinCorrelationBox, *sedPercentageBox, *bondBox, *scaleVectorBox;
    QSpinBox *atom1Box, *atom2Box, *electron1Box, *electron2Box, *eigenvectorSpinBox;
    QPushButton *deselectAllButton;
    QTreeWidget *maximaList;
    QLabel *probabilitySum, *eigenvalueLabel;
    std::vector<int> lastMovedElectronClusterVector;

    void showSplashScreen();
    void loadData();
    void initialView();
    void setupSpinBoxes();
    void setupLabels();
    void connectSignals();
    void createWidget();
    void redrawSpinDecorations();
    double sumProbabilities();
    void resetEigenvalueLabel();
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
