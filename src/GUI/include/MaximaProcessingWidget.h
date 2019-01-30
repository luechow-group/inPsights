//
// Created by heuer on 05.12.18.
//

#ifndef INPSIGHTS_MAXIMAPROCESSINGWIDGET_H
#define INPSIGHTS_MAXIMAPROCESSINGWIDGET_H

#include <QWidget>
#include <QTreeWidget>
#include <Statistics.h>
#include <QLabel>
#include <QGridLayout>

class ClusterData;

class MaximaProcessingWidget : public QWidget {
Q_OBJECT
public:
    explicit MaximaProcessingWidget(QWidget* parent = nullptr);

    void initializeTreeItems(QTreeWidget &tree, int numberOfParticles);
    void setAtomEnergies(SingleParticlesStatistics EnStats);
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
    bool initializedQ_;
    SingleParticlesStatistics EnStats_;
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
