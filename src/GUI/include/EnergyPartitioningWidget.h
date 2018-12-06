//
// Created by heuer on 05.12.18.
//

#ifndef INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
#define INPSIGHTS_ENERGYPARTITIONINGWIDGET_H

#include <QWidget>
#include <QTableWidget>

#include <Statistics.h>

#include <ClusterData.h>

class EnergyPartitioningWidget : public QWidget {
Q_OBJECT
public:
    explicit EnergyPartitioningWidget(QWidget* parent = nullptr, int nAtoms = 0, int nElectrons = 0);

    void initializeItems(int nAtoms, int nElectrons);
    void setAtomEnergies(IntraParticlesStatistics VnnStats);
    void updateData(ClusterData& clusterData) const ;

private:
    bool initializedQ_;
    IntraParticlesStatistics VnnStats_;
    QTableWidget *Ee, *En;
};


class EnergyPartitioningWidget2 : public QWidget {
Q_OBJECT
public:
    // Improve constructor
    EnergyPartitioningWidget2(const IntraParticlesStatistics& VnnStats, int nAtoms = 0, int nElectrons = 0, QWidget* parent = nullptr);

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
