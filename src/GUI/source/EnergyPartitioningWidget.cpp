#include <utility>

//
// Created by heuer on 05.12.18.
//

#include <EnergyPartitioningWidget.h>
#include <QHeaderView>
#include <ClusterData.h>
#include <OneParticleEnergies.h>
#include <IntegerSortedTreeWidgetItem.h>

EnergyPartitioningWidget::EnergyPartitioningWidget(QWidget *parent)
        :
        QWidget(parent),
        console(Logger::get()),
        initializedQ_(false),
        grid_(this) {

    grid_.addWidget(&Ee_, 0, 0, Qt::AlignCenter);
    grid_.addWidget(&En_, 0, 1, Qt::AlignCenter);
    grid_.addWidget(&Eintra_, 1, 0, Qt::AlignRight);
    grid_.addWidget(&Einter_, 2, 0, Qt::AlignRight);
    grid_.addWidget(&EintraErr_, 1, 1, Qt::AlignLeft);
    grid_.addWidget(&EinterErr_, 2, 1, Qt::AlignLeft);

    initializeTree(Ee_,QString("e"));
    initializeTree(En_,QString("n"));

    connect(&Ee_, SIGNAL(itemChanged(QTreeWidgetItem*, int)),
            this, SLOT(onItemChanged()));
    connect(&En_, SIGNAL(itemChanged(QTreeWidgetItem*, int)),
            this, SLOT(onItemChanged()));

    connect(&En_, SIGNAL(itemChanged(QTreeWidgetItem*,int)),
            this, SLOT(onAtomSelectionChanged()));
    connect(&Ee_, SIGNAL(itemChanged(QTreeWidgetItem*,int)),
            this, SLOT(onElectronSelectionChanged()));

}

void EnergyPartitioningWidget::initializeTree(QTreeWidget &tree, const QString& particleSymbol) const {
    tree.setColumnCount(3);
    tree.setHeaderLabels(QStringList({
        QString("E%1 / Eh").arg(particleSymbol),
        QString("u_E%1 / Eh").arg(particleSymbol),
        QString("ID")}));
    tree.setSortingEnabled(true);
    tree.header()->setStretchLastSection(false);
    //tree.setMinimumWidth(150);
    //tree.setSizePolicy(QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Minimum);
    tree.setFixedWidth(175);
}

QTreeWidget& EnergyPartitioningWidget::atomsTreeWidget() {
    return En_;
}

QTreeWidget& EnergyPartitioningWidget::electronsTreeWidget() {
    return Ee_;
}

void EnergyPartitioningWidget::onItemChanged() {
    recalculateMotifEnergy();
}

void EnergyPartitioningWidget::recalculateMotifEnergy(){
    double intra = 0, inter = 0, intraErr = 0, interErr = 0;
    addContributions(Ee_, intra, inter, intraErr, interErr);
    addContributions(En_, intra, inter, intraErr, interErr);

    Eintra_.setText(QString("E_intra = ") + QString::number(intra, 'f', 4));
    Einter_.setText(QString("E_inter = ") + QString::number(inter, 'f', 4));
    EintraErr_.setText(QString(" ± ") + QString::number(sqrt(intraErr), 'f', 4) + QString(" Eh"));
    EinterErr_.setText(QString(" ± ") + QString::number(sqrt(interErr), 'f', 4) + QString(" Eh"));
}


void EnergyPartitioningWidget::onAtomSelectionChanged() {
    emit atomsChecked(getSelectedItems(En_));
}

void EnergyPartitioningWidget::onElectronSelectionChanged() {
    emit electronsChecked(getSelectedItems(Ee_));
}

std::vector<int> EnergyPartitioningWidget::getSelectedItems(const QTreeWidget &tree) {
    std::vector<int> selection;

    for (int i = 0; i < tree.topLevelItemCount(); ++i) {
        auto item = tree.topLevelItem(i);

        if(item->checkState(0) == Qt::Checked) {
            auto id = item->data(2, Qt::UserRole).toInt();
            selection.push_back(id);
        }
    }


    return selection;
}

void EnergyPartitioningWidget::addContributions(const QTreeWidget &tree,
                                                double &intra, double &inter, double &intraErr, double &interErr) const {
    for (int i = 0; i < tree.topLevelItemCount(); ++i) {
        auto topItem = tree.topLevelItem(i);
        auto value = topItem->data(0, Qt::UserRole).toDouble();
        auto squaredError = std::pow(topItem->data(1, Qt::UserRole).toDouble(),2);

        if(topItem->checkState(0) == Qt::Checked) {
            intra += value;
            intraErr += squaredError;
        } else {
            inter += value;
            interErr += squaredError;
        }
    }
}

void EnergyPartitioningWidget::initializeTreeItems(QTreeWidget &tree, int numberOfParticles) {
    for (int i = 0; i < numberOfParticles; ++i) {
        auto item = new IntegerSortedTreeWidgetItem();
        item->setCheckState(0,Qt::CheckState::Unchecked);
        tree.addTopLevelItem(item);
    }
}

void EnergyPartitioningWidget::setAtomEnergies(IntraParticlesStatistics VnnStats) {
    VnnStats_ = std::move(VnnStats);
}

void EnergyPartitioningWidget::updateData(const ClusterData &clusterData) {
    assert(Ee_.topLevelItemCount() == static_cast<int>(clusterData.VenStats_.rows())
    && "The number of tree items and electrons must match.");
    assert(En_.topLevelItemCount() == static_cast<int>(clusterData.VenStats_.cols())
    && "The number of tree items and atoms must match.");

    updateEnergies(Ee_,
            OneParticleEnergies::oneElectronEnergies(clusterData),
            OneParticleEnergies::oneElectronEnergiesErrors(clusterData));

    updateEnergies(En_,
                   OneParticleEnergies::oneAtomEnergies(VnnStats_, clusterData),
                   OneParticleEnergies::oneAtomEnergiesErrors(VnnStats_, clusterData));

}

void EnergyPartitioningWidget::updateEnergies(QTreeWidget &tree,
                                              const Eigen::VectorXd &energies,
                                              const Eigen::VectorXd &errors) const {
    assert(energies.size() == errors.size()
    && "Value and error vector must have the same size.");
    assert(tree.topLevelItemCount() == static_cast<int>(energies.size())
    && "The number of tree items and energy values must match.");

    tree.setSortingEnabled(false);
    for (int i = 0; i < energies.size(); ++i) {
        auto item = tree.topLevelItem(i);
        item->setData(0, Qt::UserRole, energies[i]);
        item->setData(1, Qt::UserRole, errors[i]);
        item->setData(2, Qt::UserRole, i);
        item->setText(0, QString::number(energies[i], 'f', 4));
        item->setText(1, QString::number(errors[i], 'f', 4));
        item->setText(2, QString::number(i));
    }
    tree.setSortingEnabled(true);
    tree.resizeColumnToContents(0);
    tree.resizeColumnToContents(1);
    tree.resizeColumnToContents(2);
}

