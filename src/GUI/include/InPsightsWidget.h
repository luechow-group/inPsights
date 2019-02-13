//
// Created by Michael Heuer on 18.11.18.
//

#ifndef INPSIGHTS_INPSIGHTSWIDGET_H
#define INPSIGHTS_INPSIGHTSWIDGET_H

#include <MoleculeWidget.h>
#include <MaximaProcessingWidget.h>

#include <QCheckBox>
#include <QSlider>
#include <QTreeWidget>
#include <ClusterData.h>

class InPsightsWidget : public QWidget {
Q_OBJECT
public:
    explicit InPsightsWidget(QWidget *parent = nullptr);

public slots:
    void selectedStructure(QTreeWidgetItem *item, int column);
    void onAtomsChecked(int stateId);
    void onBondsChecked(int stateId);
    void onSpinConnectionsChecked(int stateId = 0);
    void onSpinCorrelationsChecked(int stateId= 0);
    void onSpinCorrelationsSliderChanged(int value);

private:
    MoleculeWidget *moleculeWidget;
    MaximaProcessingWidget *maximaProcessingWidget;
    QCheckBox *atomsCheckBox, *bondsCheckBox, *spinConnectionsCheckBox, *spinCorrelationsCheckBox;
    QSlider *spinCorrelationSlider;
    QLabel *spinCorrelationSliderLabel;
    QTreeWidget *maximaList;
    std::vector<ClusterData> clusterCollection_;

    void showSplashScreen();
    void loadData();
    void initialView();
    void setupSliderBox();
    void updateSpinCorrelationSliderLabel(int value);
    void connectSignals();
    void createWidget();
    void redrawSpinDecorations();
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
