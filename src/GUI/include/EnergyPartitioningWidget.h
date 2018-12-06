//
// Created by heuer on 05.12.18.
//

#ifndef INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
#define INPSIGHTS_ENERGYPARTITIONINGWIDGET_H

#include <QWidget>
#include <QTableWidget>
#include <QGridLayout>
#include <Statistics.h>

class ClusterData;

class EnergyPartitioningWidget : public QWidget {
Q_OBJECT
public:
    // Improve constructor
    EnergyPartitioningWidget(const IntraParticlesStatistics& VnnStats, int nAtoms = 0, int nElectrons = 0, QWidget* parent = nullptr);

    void setAtomEnergies(IntraParticlesStatistics VnnStats);
    void updateData(ClusterData& clusterData) const ;
    void setTableSizes(int nAtoms = 0, int nElectrons = 0) const;

private:
    IntraParticlesStatistics VnnStats_;
    QTableWidget *Te, *Ven, *Vee, *Ee, *Vnn, *En;
    QList<QTableWidget*> tables;

    void adjustAllTableSizes() const;

    void placeItem(QTableWidget *table, double value, int i, int j = 0) const;

    void setTableSize(QTableWidget *table, int rows, int cols) const;
};

#endif //INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
