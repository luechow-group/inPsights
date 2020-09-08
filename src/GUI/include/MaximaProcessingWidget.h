// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MAXIMAPROCESSINGWIDGET_H
#define INPSIGHTS_MAXIMAPROCESSINGWIDGET_H

#include <QWidget>
#include <QTreeWidget>
#include <Statistics.h>
#include <QLabel>
#include <QGridLayout>
#include <ParticlesVector.h>

class ClusterData;

class MaximaProcessingWidget : public QWidget {
Q_OBJECT
public:
    explicit MaximaProcessingWidget(QWidget* parent = nullptr);

    void initializeTreeItems(QTreeWidget &tree, int numberOfParticles);
    void setAtomsVector(const AtomsVector& atoms);
    void setAtomEnergies(VectorStatistics EnStats);
    void updateData(const ClusterData& clusterData);

    QTreeWidget& atomsTreeWidget();
    QTreeWidget& electronsTreeWidget();

Q_SIGNALS:
    void atomsChecked(std::vector<int> selectedIds);
    void electronsChecked(std::vector<int> selectedIds);
    void atomsHighlighted(std::vector<int> selectedIds);
    void electronsHighlighted(std::vector<int> selectedIds);

public Q_SLOTS:
    void onAtomItemChanged();
    void onElectronItemChanged();
    void onAtomSelectionChanged();
    void onElectronSelectionChanged();

private:
    AtomsVector atoms_;
    bool initializedQ_;
    VectorStatistics EnStats_;
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
    std::vector<int> getCheckedItems(const QTreeWidget &tree);
    std::vector<int> getSelectedItems(const QTreeWidget &tree);
};

#endif //INPSIGHTS_MAXIMAPROCESSINGWIDGET_H
