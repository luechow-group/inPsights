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
    void onSpinConnectionsChecked(int stateId = 0);
    void onSpinCorrelationsChecked(int stateId= 0);
    void onPlotAllChecked(int stateId= 0);
    void onSpinCorrelationsBoxChanged(double value);

private:
    std::string filename_;
    MoleculeWidget *moleculeWidget;
    MaximaProcessingWidget *maximaProcessingWidget;
    QCheckBox *atomsCheckBox, *bondsCheckBox, *axesCheckBox, *sampleAverageCheckBox, *spinConnectionsCheckBox, *spinCorrelationsCheckBox,
    *sedsCheckBox,*maximaHullsCheckBox, *plotAllCheckBox, *coloredCheckBox;
    QDoubleSpinBox *spinCorrelationBox, *sedPercentageBox;
    QTreeWidget *maximaList;
    std::vector<ClusterData> clusterCollection_;

    void showSplashScreen();
    void loadData();
    void initialView();
    void setupSpinBoxes();
    void connectSignals();
    void createWidget();
    void redrawSpinDecorations();
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
