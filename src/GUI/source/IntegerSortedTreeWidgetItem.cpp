//
// Created by heuer on 05.12.18.
//

#include <IntegerSortedTreeWidgetItem.h>

#include "IntegerSortedTreeWidgetItem.h"

IntegerSortedTreeWidgetItem::IntegerSortedTreeWidgetItem()
        : QTreeWidgetItem() {}

IntegerSortedTreeWidgetItem::IntegerSortedTreeWidgetItem(QTreeWidget *tree)
        : QTreeWidgetItem(tree) {}

IntegerSortedTreeWidgetItem::IntegerSortedTreeWidgetItem(QTreeWidget *parent, const QStringList &strings)
        : QTreeWidgetItem(parent, strings) {}

IntegerSortedTreeWidgetItem::IntegerSortedTreeWidgetItem(QTreeWidgetItem *parent, const QStringList &strings, int type)
        : QTreeWidgetItem(parent, strings, type) {}

bool IntegerSortedTreeWidgetItem::operator<(const QTreeWidgetItem &other) const {
    auto sortCol = treeWidget()->sortColumn();
    auto myNumber = text(sortCol).toDouble(); //TODO Refactor
    auto  otherNumber = other.text(sortCol).toDouble();
    return myNumber < otherNumber;
}