//
// Created by Michael Heuer on 18.11.18.
//

#ifndef INPSIGHTS_INPSIGHTSWIDGET_H
#define INPSIGHTS_INPSIGHTSWIDGET_H

#include <Logger.h>
#include <MoleculeWidget.h>
#include <QWidget>

#include <QCheckBox>
#include <QSlider>
#include <QTreeWidget>
#include <ClusterData.h>

class QSplashScreen;

class InPsightsWidget : public QWidget {
Q_OBJECT
public:
    explicit InPsightsWidget(QWidget *parent = nullptr);

public slots:

    void selectedStructure(QTreeWidgetItem *item, int column);

    void onAtomsChecked(int stateId);

    void onBondsChecked(int stateId);

    void onSpinConnectionsChecked(int stateId);

    void onSpinCorrelationsChecked(int stateId);

    void onSpinCorrelationsSliderChanged(int value);

private:
    std::shared_ptr<spdlog::logger> console;
    MoleculeWidget *moleculeWidget;
    QCheckBox *atomsCheckBox, *bondsCheckBox, *spinConnectionsCheckBox, *spinCorrelationsCheckBox;
    QSlider *spinCorrelationSlider;
    QLabel *spinCorrelationSliderLabel;
    QTreeWidget *maximaList;
    std::vector<ClusterData> clusterCollection_;

    QSplashScreen *createSplashScreen();

    void loadData();

    void initialView();

    void setupSliderBox();

    void updateSpinCorrelationSliderLabel(int value);
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
