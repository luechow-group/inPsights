//
// Created by heuer on 05.12.18.
//

#ifndef INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
#define INPSIGHTS_ENERGYPARTITIONINGWIDGET_H



#include <QWidget>
#include <QTreeWidget>
#include <Statistics.h>
#include <QLabel>
#include <QGridLayout>
#include <Logger.h>

class ClusterData;

class EnergyPartitioningWidget : public QWidget {
Q_OBJECT
public:
    explicit EnergyPartitioningWidget(QWidget* parent = nullptr);

    void initializeTreeItems(QTreeWidget &tree, int numberOfParticles);
    void setAtomEnergies(IntraParticlesStatistics VnnStats);
    void updateData(const ClusterData& clusterData);

    QTreeWidget& atomsTreeWidget();
    QTreeWidget& electronsTreeWidget();

Q_SIGNALS:
    void atomsChecked(std::vector<int> selectedIds);
    void electronsChecked(std::vector<int> selectedIds);

public Q_SLOTS:
    void onItemChanged();
    void onAtomSelectionChanged();
    void onElectronSelectionChanged();

private:
    std::shared_ptr<spdlog::logger> console;
    bool initializedQ_;
    IntraParticlesStatistics VnnStats_;
    QGridLayout grid_;
    QTreeWidget Ee_, En_;
    QLabel Eintra_, Einter_,EintraErr_, EinterErr_;

    void updateEnergies(QTreeWidget &tree,
                        const Eigen::VectorXd &energies,
                        const Eigen::VectorXd &errors) const;

    void initializeTree(QTreeWidget &tree, const QString& particleSymbol) const;
    //TODO use struct
    void addContributions(const QTreeWidget &tree, double &intra, double &inter, double &intraErr, double &interErr) const;
    void recalculateMotifEnergy();
    std::vector<int> getSelectedItems(const QTreeWidget &tree);
};

#endif //INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
