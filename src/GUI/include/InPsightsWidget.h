//
// Created by Michael Heuer on 18.11.18.
//

#ifndef INPSIGHTS_INPSIGHTSWIDGET_H
#define INPSIGHTS_INPSIGHTSWIDGET_H

#include <MoleculeWidget.h>
#include <QWidget>

#include <QCheckBox>
#include <QSlider>
#include <QListWidget>
#include <ClusterData.h>

class QSplashScreen;

class InPsightsWidget : public QWidget {
Q_OBJECT
public:
    explicit InPsightsWidget(QWidget *parent = nullptr);

public slots:
    void selectedStructure(QListWidgetItem* item);

    void onBondsChecked();

    void onConnectionsChecked();

private:
    MoleculeWidget *moleculeWidget_;
    QCheckBox *bondsCheckBox_, *spinConnectionsCheckBox_, *spinCorrelationsCheckBox_;
    QSlider *spinCorrelationSlider_;
    QLabel *spinCorrelationSliderLabel_;
    QListWidget *maximaList_;
    std::vector<ClusterData> clusterCollection_;

    QSplashScreen *createSplashScreen();

    void loadData();

    void initialView();
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
