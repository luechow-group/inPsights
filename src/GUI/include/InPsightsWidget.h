//
// Created by Michael Heuer on 18.11.18.
//

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

public slots:
    void selectedStructure(QTreeWidgetItem *item, int column);
    void onAtomsChecked(int stateId);
    void onBondsChecked(int stateId);
    void onAxesChecked(int stateId);
    void onSpinConnectionsChecked(int stateId = 0);
    void onSpinCorrelationsChecked(int stateId= 0);
    void onSpinCorrelationsBoxChanged(double value);

private:
    std::string filename_;
    MoleculeWidget *moleculeWidget;
    MaximaProcessingWidget *maximaProcessingWidget;
    QCheckBox *atomsCheckBox, *bondsCheckBox, *axesCheckBox, *spinConnectionsCheckBox, *spinCorrelationsCheckBox, *sedsCheckBox;
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
