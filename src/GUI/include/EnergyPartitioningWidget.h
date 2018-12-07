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

    void initializeTreeItems(QTreeWidget *tree, int numberOfParticles);
    void setAtomEnergies(IntraParticlesStatistics VnnStats);
    void updateData(const ClusterData& clusterData);

    QTreeWidget* atomsTreeWidget();
    QTreeWidget* electronsTreeWidget();

private:
    bool initializedQ_;
    IntraParticlesStatistics VnnStats_;
    QTreeWidget *Ee, *En;

    void updateEnergies(QTreeWidget *tree,
                        const Eigen::VectorXd &energies,
                        const Eigen::VectorXd &errors) const;

    void initializeTree(QTreeWidget *tree) const;
};

#endif //INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
