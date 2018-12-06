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

EnergyPartitioningWidget::EnergyPartitioningWidget(QWidget *parent, int nAtoms, int nElectrons)
        :
        QWidget(parent),
        initializedQ_(false),
        Ee(new QTableWidget(nElectrons, 2, this)),
        En(new QTableWidget(nAtoms, 2, this)) {

    auto outerLayout = new QHBoxLayout(this);

    outerLayout->addWidget(Ee);
    outerLayout->addWidget(En);

    QFont font;
    font.setPointSize(8);
    Ee->setFont(font);
    En->setFont(font);
}

void EnergyPartitioningWidget::initializeItems(int nAtoms, int nElectrons) {
    Ee->setRowCount(nElectrons);
    Ee->setColumnCount(2);
    for (int i = 0; i < nElectrons; ++i) {
        Ee->setItem(i, 0, new QTableWidgetItem());
        Ee->setItem(i, 1, new QTableWidgetItem());
    }

    En->setRowCount(nAtoms);
    En->setColumnCount(2);
    for (int i = 0; i < nAtoms; ++i) {
        En->setItem(i, 0, new QTableWidgetItem());
        En->setItem(i, 1, new QTableWidgetItem());
    }

}

void EnergyPartitioningWidget::setAtomEnergies(IntraParticlesStatistics VnnStats) {
    VnnStats_ = std::move(VnnStats);
}

void EnergyPartitioningWidget::updateData(ClusterData &clusterData) const {
    auto nElectrons = static_cast<int>(clusterData.VenStats_.rows());
    auto nAtoms = static_cast<int>(clusterData.VenStats_.cols());

    Ee->setRowCount(nElectrons);
    En->setRowCount(nAtoms);

    auto EeVal = OneParticleEnergies::oneElectronEnergies(clusterData);
    auto EeErr = OneParticleEnergies::oneElectronEnergiesErrors(clusterData);

    auto EnVal = OneParticleEnergies::oneAtomEnergies(VnnStats_,clusterData);
    auto EnErr = OneParticleEnergies::oneAtomEnergiesErrors(VnnStats_,clusterData);

    for (int i = 0; i < nElectrons; ++i) {
        Ee->item(i, 0)->setData(Qt::UserRole, EeVal[i]);
        Ee->item(i, 1)->setData(Qt::UserRole, EeErr[i]);
        Ee->item(i, 0)->setText(QString::number(EeVal[i], 'f', 4));
        Ee->item(i, 1)->setText(QString::number(EeErr[i], 'f', 4));
    }

    for (int i = 0; i < nAtoms; ++i) {
        En->item(i, 0)->setData(Qt::UserRole, EnVal[i]);
        En->item(i, 1)->setData(Qt::UserRole, EnErr[i]);
        En->item(i, 0)->setText(QString::number(EnVal[i], 'f', 4));
        En->item(i, 1)->setText(QString::number(EnErr[i], 'f', 4));
    }
}

EnergyPartitioningWidget2::EnergyPartitioningWidget2(const IntraParticlesStatistics& VnnStats, int nAtoms, int nElectrons, QWidget *parent)
        : QWidget(parent),
          VnnStats_(VnnStats),
          Te(new QTableWidget(1, nElectrons, this)),
          Ven(new QTableWidget(nAtoms, nElectrons, this)),
          Vee(new QTableWidget(nElectrons, nElectrons, this)),
          Ee(new QTableWidget(1, nElectrons, this)),
          Vnn(new QTableWidget(nAtoms, nAtoms, this)),
          En(new QTableWidget(nAtoms, 1, this)),
          tables({Te, Ven, Vee, Ee, Vnn, En}) {

    auto gridLayout = new QGridLayout(this);

    gridLayout->addWidget(Te,  0,0);
    gridLayout->addWidget(Vee, 1,0);
    gridLayout->addWidget(Ven, 2,0);
    gridLayout->addWidget(Ee,  3,0);
    gridLayout->addWidget(Vnn, 2,1);
    gridLayout->addWidget(En,  2,2);

    setTableSizes(nAtoms,nElectrons);

    QFont font;
    font.setPointSize(8);
    for(auto & t : tables)
        t->setFont(font);
}

void EnergyPartitioningWidget2::setTableSize(QTableWidget *table, int rows, int cols) const {
    table->setRowCount(rows);
    table->setColumnCount(cols);

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            table->setItem(i, j, new QTableWidgetItem());
}

void EnergyPartitioningWidget2::setTableSizes(int nAtoms, int nElectrons) const {
    setTableSize(Te , 1         , nElectrons);
    setTableSize(Ven, nAtoms    , nElectrons);
    setTableSize(Vee, nElectrons, nElectrons);
    setTableSize(Ee , 1         , nElectrons);
    setTableSize(Vnn, nAtoms    , nAtoms    );
    setTableSize(En , nAtoms    , 1         );
}

void EnergyPartitioningWidget2::adjustAllTableSizes() const {
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


void EnergyPartitioningWidget2::updateData(ClusterData &clusterData) const {

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

void EnergyPartitioningWidget2::placeItem(QTableWidget *table, double value, int j, int i) const {
    table->item(i,j)->setData(Qt::UserRole, value);
    table->item(i,j)->setText(QString::number(value, 'f', 4));
}
