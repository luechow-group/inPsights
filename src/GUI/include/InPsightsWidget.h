//
// Created by Michael Heuer on 18.11.18.
//

#ifndef INPSIGHTS_INPSIGHTSWIDGET_H
#define INPSIGHTS_INPSIGHTSWIDGET_H

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
    void selectedStructure(QTreeWidgetItem* item, int column);

    void onAtomsChecked(int stateId);

    void onBondsChecked(int stateId);

    void onSpinConnectionsChecked(int stateId);

    void onSpinCorrelationsChecked(int stateId);

    void onSpinCorrelationsSliderChanged(int value);

private:
    MoleculeWidget *moleculeWidget_;
    QCheckBox *atomsCheckBox_, *bondsCheckBox_, *spinConnectionsCheckBox_, *spinCorrelationsCheckBox_;
    QSlider *spinCorrelationSlider_;
    QLabel *spinCorrelationSliderLabel_;
    QTreeWidget *maximaList_;
    std::vector<ClusterData> clusterCollection_;

    QSplashScreen *createSplashScreen();

    void loadData();

    void initialView();
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
