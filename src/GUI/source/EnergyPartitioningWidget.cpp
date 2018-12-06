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

    Ee->setColumnCount(2);
    Ee->setHeaderLabels(QList<QString>({"Ee", "Error"}));
    Ee->setSortingEnabled(true);
    Ee->header()->setStretchLastSection(false);
    //Ee->setMinimumWidth(150);
    //Ee->setSizePolicy(QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Minimum);
    Ee->setFixedWidth(150);

    En->setColumnCount(2);
    En->setHeaderLabels(QList<QString>({"En", "Error"}));
    En->setSortingEnabled(true);
    En->header()->setStretchLastSection(false);
    //En->setMinimumWidth(150);
    //En->setSizePolicy(QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Minimum);
    En->setFixedWidth(150);
}

void EnergyPartitioningWidget::initializeItems(int nAtoms, int nElectrons) {
    for (int i = 0; i < nElectrons; ++i)
        Ee->addTopLevelItem(new QTreeWidgetItem());

    for (int i = 0; i < nAtoms; ++i)
        En->addTopLevelItem(new QTreeWidgetItem());

}

void EnergyPartitioningWidget::setAtomEnergies(IntraParticlesStatistics VnnStats) {
    VnnStats_ = std::move(VnnStats);
}

void EnergyPartitioningWidget::updateData(const ClusterData &clusterData) {
    auto nElectrons = static_cast<int>(clusterData.VenStats_.rows());
    auto nAtoms = static_cast<int>(clusterData.VenStats_.cols());

    auto EeVal = OneParticleEnergies::oneElectronEnergies(clusterData);
    auto EeErr = OneParticleEnergies::oneElectronEnergiesErrors(clusterData);

    auto EnVal = OneParticleEnergies::oneAtomEnergies(VnnStats_,clusterData);
    auto EnErr = OneParticleEnergies::oneAtomEnergiesErrors(VnnStats_,clusterData);

    for (int i = 0; i < nElectrons; ++i) {
        auto item = Ee->topLevelItem(i);
        item->setData(0,Qt::UserRole, EeVal[i]);
        item->setData(1,Qt::UserRole, EeErr[i]);
        item->setText(0,QString::number(EeVal[i],'f',4));
        item->setText(1,QString::number(EeErr[i],'f',4));
    }

    for (int k = 0; k < nAtoms; ++k) {
        auto item = En->topLevelItem(k);
        item->setData(0,Qt::UserRole, EnVal[k]);
        item->setData(1,Qt::UserRole, EnErr[k]);
        item->setText(0,QString::number(EnVal[k],'f',4));
        item->setText(1,QString::number(EnErr[k],'f',4));
    }

    Ee->resizeColumnToContents(0);
    Ee->resizeColumnToContents(1);

    En->resizeColumnToContents(0);
    En->resizeColumnToContents(1);

    //this->resize(minimumWidth(), sizeHint().height());
}

