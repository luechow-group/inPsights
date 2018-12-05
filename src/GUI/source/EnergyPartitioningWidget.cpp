//
// Created by heuer on 05.12.18.
//

#include <EnergyPartitioningWidget.h>
#include <QHeaderView>

EnergyPartitioningWidget::EnergyPartitioningWidget(QWidget *parent, int nAtoms, int nElectrons)
        : QWidget(parent),
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

    for(auto t : tables) {
        t->setSizePolicy(QSizePolicy(QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Minimum));
        t->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        t->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        t->resizeColumnsToContents();

        //t->setFixedSize(20,40);

        //auto hh= t->horizontalHeader();
        //auto vh= t->verticalHeader()
        t->setFixedSize(
                t->horizontalHeader()->length()+ t->verticalHeader()->width(),
                t->verticalHeader()->height() + t->horizontalHeader()->height());
    }
}