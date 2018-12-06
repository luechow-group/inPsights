//
// Created by heuer on 05.12.18.
//

#ifndef INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
#define INPSIGHTS_ENERGYPARTITIONINGWIDGET_H

#include <QWidget>
#include <QTableWidget>
#include <QTreeWidget>
#include <Statistics.h>

#include <ClusterData.h>

class EnergyPartitioningWidget : public QWidget {
Q_OBJECT
public:
    explicit EnergyPartitioningWidget(QWidget* parent = nullptr);

    void initializeItems(int nAtoms, int nElectrons);
    void setAtomEnergies(IntraParticlesStatistics VnnStats);
    void updateData(const ClusterData& clusterData);

private:
    bool initializedQ_;
    IntraParticlesStatistics VnnStats_;
    QTreeWidget *Ee, *En;
};

#endif //INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
