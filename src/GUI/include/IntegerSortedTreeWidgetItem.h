//
// Created by heuer on 05.12.18.
//

#ifndef INPSIGHTS_INTEGERSORTEDTREEWIDGETITEM_H
#define INPSIGHTS_INTEGERSORTEDTREEWIDGETITEM_H

#include <QTreeWidget>

class IntegerSortedTreeWidgetItem : public QTreeWidgetItem {
public:
    IntegerSortedTreeWidgetItem();

    explicit IntegerSortedTreeWidgetItem(QTreeWidget *tree);

    IntegerSortedTreeWidgetItem(QTreeWidgetItem *parent, const QStringList &strings, int type = Type);

    IntegerSortedTreeWidgetItem(QTreeWidget * parent, const QStringList & strings);

    bool operator< (const QTreeWidgetItem &other) const override;
};

#endif //INPSIGHTS_INTEGERSORTEDTREEWIDGETITEM_H
