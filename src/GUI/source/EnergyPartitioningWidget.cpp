#include <utility>

//
// Created by heuer on 05.12.18.
//

#include <EnergyPartitioningWidget.h>
#include <QHeaderView>
#include <ClusterData.h>
#include <OneParticleEnergies.h>

EnergyPartitioningWidget::EnergyPartitioningWidget(const IntraParticlesStatistics& VnnStats, int nAtoms, int nElectrons, QWidget *parent)
        : QWidget(parent),
          VnnStats_(VnnStats),
          Te(new QTableWidget(nElectrons, 1, this)),
          Ven(new QTableWidget(nElectrons, nAtoms, this)),
          Vee(new QTableWidget(nElectrons, nElectrons, this)),
          Ee(new QTableWidget(nElectrons, 1, this)),
          Vnn(new QTableWidget(nAtoms, nAtoms, this)),
          En(new QTableWidget(1,nAtoms, this)),
          tables({Te, Ven, Vee, Ee, Vnn, En}) {

    auto gridLayout = new QGridLayout(this);

    gridLayout->addWidget(Te, 0, 0);
    gridLayout->addWidget(Ven, 0, 1);
    gridLayout->addWidget(Vee, 0, 2);
    gridLayout->addWidget(Ee, 0, 3);
    gridLayout->addWidget(Vnn, 1, 1);
    gridLayout->addWidget(En, 2, 1);

    setTableSizes(nAtoms,nElectrons);

    QFont font;
    font.setPointSize(8);
    for(auto & t : tables)
        t->setFont(font);
}

void EnergyPartitioningWidget::setTableSize(QTableWidget *table, int rows, int cols) const {
    table->setRowCount(rows);
    table->setColumnCount(cols);

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            table->setItem(i, j, new QTableWidgetItem());
}

void EnergyPartitioningWidget::setTableSizes(int nAtoms, int nElectrons) const {
    setTableSize(Te, nElectrons, 1);
    setTableSize(Ven, nElectrons, nAtoms);
    setTableSize(Vee, nElectrons, nElectrons);
    setTableSize(Ee, nElectrons, 1);
    setTableSize(Vnn, nAtoms, nAtoms);
    setTableSize(En, 1, nAtoms);
}

void EnergyPartitioningWidget::adjustAllTableSizes() const {
    for(auto t : tables) {
        t->setSizePolicy(QSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum));
        t->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        t->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        t->resizeColumnsToContents();

        t->setFixedSize(
                t->horizontalHeader()->length()+ t->verticalHeader()->width(),
                t->verticalHeader()->height() + t->horizontalHeader()->height());
    }
}


void EnergyPartitioningWidget::updateData(ClusterData &clusterData) const {

    auto nElectrons = Vee->rowCount();
    auto nAtoms = Vnn->rowCount();

    auto EeVal = OneParticleEnergies::oneElectronEnergies(clusterData);

    auto EnVal = OneParticleEnergies::oneAtomEnergies(VnnStats_, clusterData);

    for (int i = 0; i < nElectrons; ++i) {
        placeItem(Te,clusterData.TeStats_.mean()(i,0), i);
        placeItem(Ee, EeVal[i], i);

        for (int k = 0; k < nAtoms; ++k)
            placeItem(Ven, clusterData.VenStats_.mean()(i,k), i,k);


        for (int j = i+1; j < nElectrons; ++j)
            placeItem(Vee, clusterData.VeeStats_.mean()(i,j), i,j);
    }

    for (int k = 0; k < nAtoms; ++k) {
        placeItem(En, EnVal(k,0), 0, k);
        for (int l = k + 1; l < nAtoms; ++l)
            placeItem(Vnn, VnnStats_.mean()(k, l), k, l);
    }

    adjustAllTableSizes();
}

void EnergyPartitioningWidget::placeItem(QTableWidget *table, double value, int i, int j) const {
    table->item(i,j)->setData(Qt::UserRole, value);
    table->item(i,j)->setText(QString::number(value, 'f', 4));
}
