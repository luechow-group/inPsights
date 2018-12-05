//
// Created by heuer on 05.12.18.
//

#ifndef INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
#define INPSIGHTS_ENERGYPARTITIONINGWIDGET_H

#include <QWidget>
#include <QTableWidget>
#include <QGridLayout>

class EnergyPartitioningWidget : public QWidget {
Q_OBJECT
public:
    explicit EnergyPartitioningWidget(QWidget* parent = nullptr, int nAtoms = 3, int nElectrons = 5);

    QTableWidget *Te, *Ven, *Vee, *Ee, *Vnn, *En;
    QList<QTableWidget*> tables;
};

#endif //INPSIGHTS_ENERGYPARTITIONINGWIDGET_H
