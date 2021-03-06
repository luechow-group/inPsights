// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <utility>
#include <MaximaProcessingWidget.h>
#include <QHeaderView>
#include <ClusterData.h>
#include <IntegerSortedTreeWidgetItem.h>
#include <EnergyPartitioning.h>
#include <CoulombPotential.h>

MaximaProcessingWidget::MaximaProcessingWidget(QWidget *parent)
        :
        QWidget(parent),
        initializedQ_(false),
        grid_(this) {

    grid_.addWidget(&Ee_, 0, 0, Qt::AlignCenter);
    grid_.addWidget(&En_, 0, 1, Qt::AlignCenter);
    grid_.addWidget(&Eintra_, 1, 0, Qt::AlignRight);
    grid_.addWidget(&Einter_, 2, 0, Qt::AlignRight);
    grid_.addWidget(&EintraErr_, 1, 1, Qt::AlignLeft);
    grid_.addWidget(&EinterErr_, 2, 1, Qt::AlignLeft);


    initializeTree(En_,QString("n"));
    connect(&En_, &QTreeWidget::itemSelectionChanged,
            this, &MaximaProcessingWidget::onAtomSelectionChanged);
    connect(&En_, &QTreeWidget::itemChanged,
            this, &MaximaProcessingWidget::onAtomItemChanged);

    initializeTree(Ee_,QString("e"));
    connect(&Ee_, &QTreeWidget::itemSelectionChanged,
            this, &MaximaProcessingWidget::onElectronSelectionChanged);
    connect(&Ee_, SIGNAL(itemChanged(QTreeWidgetItem*,int)),
            this, SLOT(onElectronItemChanged()));
}

void MaximaProcessingWidget::setAtomsVector(const AtomsVector &atoms) {
    atoms_ = atoms;
}

void MaximaProcessingWidget::initializeTree(QTreeWidget &tree, const QString& particleSymbol) const {
    tree.setColumnCount(3);
    tree.setHeaderLabels(QStringList({
        QString("E%1 / Eh").arg(particleSymbol),
        QString("u_E%1 / Eh").arg(particleSymbol),
        QString("ID")}));
    tree.setSortingEnabled(true);
    tree.sortItems(0,Qt::SortOrder::AscendingOrder);
    tree.header()->setStretchLastSection(false);
}

QTreeWidget& MaximaProcessingWidget::atomsTreeWidget() {
    return En_;
}

QTreeWidget& MaximaProcessingWidget::electronsTreeWidget() {
    return Ee_;
}

void MaximaProcessingWidget::recalculateMotifEnergy(){
    double intra = 0, inter = 0, intraErr = 0, interErr = 0;
    addContributions(Ee_, intra, inter, intraErr, interErr);
    addContributions(En_, intra, inter, intraErr, interErr);

    Eintra_.setText(QString("E_selected = ") + QString::number(intra, 'f', 4));
    Einter_.setText(QString("E_unselected = ") + QString::number(inter, 'f', 4));
    EintraErr_.setText(QString(" ± ") + QString::number(sqrt(intraErr), 'f', 4) + QString(" Eh"));
    EinterErr_.setText(QString(" ± ") + QString::number(sqrt(interErr), 'f', 4) + QString(" Eh"));
}


void MaximaProcessingWidget::onAtomItemChanged() {
    recalculateMotifEnergy();
    emit atomsChecked(getCheckedItems(En_));
}

void MaximaProcessingWidget::onElectronItemChanged() {
    recalculateMotifEnergy();
    emit electronsChecked(getCheckedItems(Ee_));
}

void MaximaProcessingWidget::onAtomSelectionChanged() {
    emit atomsHighlighted(getSelectedItems(En_));
}

void MaximaProcessingWidget::onElectronSelectionChanged() {
    emit electronsHighlighted(getSelectedItems(Ee_));
}

std::vector<int> MaximaProcessingWidget::getSelectedItems(const QTreeWidget &tree) {
    std::vector<int> selection;

    for(auto &item : tree.selectedItems()) {
        auto id = item->data(2, Qt::UserRole).toInt();
        selection.push_back(id);
    }
    return selection;
}

std::vector<int> MaximaProcessingWidget::getCheckedItems(const QTreeWidget &tree) {
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

void MaximaProcessingWidget::addContributions(const QTreeWidget &tree,
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

void MaximaProcessingWidget::initializeTreeItems(QTreeWidget &tree, int numberOfParticles) {
    for (int i = 0; i < numberOfParticles; ++i) {
        auto item = new IntegerSortedTreeWidgetItem();
        item->setCheckState(0,Qt::CheckState::Unchecked);
        tree.addTopLevelItem(item);
    }
}

void MaximaProcessingWidget::setAtomEnergies(VectorStatistics EnStats) {
    EnStats_ = std::move(EnStats);
}


void MaximaProcessingWidget::updateData(const ClusterData &clusterData) {
    assert(Ee_.topLevelItemCount() == static_cast<int>(clusterData.EeStats_.rows())
    && "The number of tree items and electrons must match.");
    assert(En_.topLevelItemCount() == static_cast<int>(clusterData.electronicEnergyStats_.Ven().cols())
    && "The number of tree items and atoms must match.");

    auto adjacencyMatrix = GraphAnalysis::filter(clusterData.SeeStats_.mean().cwiseAbs());
    auto electrons = clusterData.representativeStructure();
    auto Vnn = CoulombPotential::energies(atoms_);

    updateEnergies(Ee_, clusterData.EeStats_.mean(), clusterData.EeStats_.standardError());
    updateEnergies(En_, EnStats_.mean(), EnStats_.standardError());
}

void MaximaProcessingWidget::updateEnergies(QTreeWidget &tree,
                                              const Eigen::VectorXd &energies,
                                              const Eigen::VectorXd &errors) const {
    assert(energies.size() == errors.size()
    && "Value and error vector must have the same size.");
    assert(tree.topLevelItemCount() == static_cast<int>(energies.size())
    && "The number of tree items and energy values must match.");

    tree.setSortingEnabled(false);
    for (Eigen::Index i = 0; i < energies.size(); ++i) {
        auto item = tree.topLevelItem(i);
        item->setData(0, Qt::UserRole, energies[i]);
        item->setData(1, Qt::UserRole, errors[i]);
        item->setData(2, Qt::UserRole, int(i));
        item->setText(0, QString::number(energies[i], 'f', 4));
        item->setText(1, QString::number(errors[i], 'f', 4));
        item->setText(2, QString::number(i));
    }
    tree.setSortingEnabled(true);
    tree.resizeColumnToContents(0);
    tree.resizeColumnToContents(1);
    tree.resizeColumnToContents(2);
}
