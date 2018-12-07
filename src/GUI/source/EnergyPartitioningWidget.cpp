#include <utility>

//
// Created by heuer on 05.12.18.
//

#include <EnergyPartitioningWidget.h>
#include <QHeaderView>
#include <ClusterData.h>
#include <OneParticleEnergies.h>
#include <QGridLayout>
#include <QHBoxLayout>

EnergyPartitioningWidget::EnergyPartitioningWidget(QWidget *parent)
        :
        QWidget(parent),
        initializedQ_(false),
        Ee(new QTreeWidget(this)),
        En(new QTreeWidget(this)) {

    auto outerLayout = new QHBoxLayout(this);

    outerLayout->addWidget(Ee);
    outerLayout->addWidget(En);

    initializeTree(Ee);
    initializeTree(En);
}

void EnergyPartitioningWidget::initializeTree(QTreeWidget *tree) const {
    tree->setColumnCount(3);
    tree->setHeaderLabels(QList<QString>({"Ee", "Error","ID"}));
    tree->setSortingEnabled(true);
    tree->header()->setStretchLastSection(false);
    //tree->setMinimumWidth(150);
    //tree->setSizePolicy(QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Minimum);
    tree->setFixedWidth(150);
}

QTreeWidget *EnergyPartitioningWidget::atomsTreeWidget() {
    return En;
}

QTreeWidget *EnergyPartitioningWidget::electronsTreeWidget() {
    return Ee;

}

void EnergyPartitioningWidget::initializeTreeItems(QTreeWidget *tree, int numberOfParticles) {
    for (int i = 0; i < numberOfParticles; ++i) {
        auto item = new QTreeWidgetItem();
        item->setCheckState(0,Qt::CheckState::Unchecked);
        tree->addTopLevelItem(item);
    }
}

void EnergyPartitioningWidget::setAtomEnergies(IntraParticlesStatistics VnnStats) {
    VnnStats_ = std::move(VnnStats);
}

void EnergyPartitioningWidget::updateData(const ClusterData &clusterData) {
    assert(Ee->topLevelItemCount() == static_cast<int>(clusterData.VenStats_.rows())
    && "The number of tree items and electrons must match.");
    assert(En->topLevelItemCount() == static_cast<int>(clusterData.VenStats_.cols())
    && "The number of tree items and atoms must match.");

    updateEnergies(Ee,
            OneParticleEnergies::oneElectronEnergies(clusterData),
            OneParticleEnergies::oneElectronEnergiesErrors(clusterData));

    updateEnergies(En,
                   OneParticleEnergies::oneAtomEnergies(VnnStats_, clusterData),
                   OneParticleEnergies::oneAtomEnergiesErrors(VnnStats_, clusterData));

}

void EnergyPartitioningWidget::updateEnergies(QTreeWidget *tree,
                                              const Eigen::VectorXd &energies,
                                              const Eigen::VectorXd &errors) const {
    assert(energies.size() == errors.size()
    && "Value and error vector must have the same size.");
    assert(tree->topLevelItemCount() == static_cast<int>(energies.size())
    && "The number of tree items and energy values must match.");

    for (int i = 0; i < energies.size(); ++i) {
        auto item = tree->topLevelItem(i);
        item->setData(0, Qt::UserRole, energies[i]);
        item->setData(1, Qt::UserRole, errors[i]);
        item->setData(2, Qt::UserRole, i);
        item->setText(0, QString::number(energies[i], 'f', 4));
        item->setText(1, QString::number(errors[i], 'f', 4));
        item->setText(2, QString::number(i));
    }
    tree->resizeColumnToContents(0);
    tree->resizeColumnToContents(1);
    tree->resizeColumnToContents(2);
}

