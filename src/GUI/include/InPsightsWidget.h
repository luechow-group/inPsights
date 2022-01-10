// Copyright (C) 2018-2019 Michael Heuer.
// Copyright (C) 2021 Leonard Reuter.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_INPSIGHTSWIDGET_H
#define INPSIGHTS_INPSIGHTSWIDGET_H

#include <MoleculeWidget.h>
#include <MaximaProcessingWidget.h>

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
    void onSedsExport(bool);
    void onDeselectAll(bool);
    void onSedChecked(int stateId);
    void updateSelectedStructures(int stateId);

private:
    std::string filename_;
    MoleculeWidget *moleculeWidget;
    //MaximaProcessingWidget *maximaProcessingWidget;
    QCheckBox *atomsCheckBox, *bondsCheckBox, *axesCheckBox, *sampleAverageCheckBox, *spinCorrelationsCheckBox,
    *sedsCheckBox,*maximaHullsCheckBox, *plotAllCheckBox, *coloredCheckBox;
    QDoubleSpinBox *spinCorrelationBox, *sedPercentageBox;
    QSpinBox *atom1Box, *atom2Box, *electron1Box, *electron2Box;
    QPushButton *sedsExportButton, *deselectAllButton;
    QTreeWidget *maximaList;
    QLabel *probabilitySum;
    std::vector<ClusterData> clusterCollection_;

    void showSplashScreen();
    void loadData();
    void initialView();
    void setupSpinBoxes();
    void connectSignals();
    void createWidget();
    void redrawSpinDecorations();
    double sumProbabilities();
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
