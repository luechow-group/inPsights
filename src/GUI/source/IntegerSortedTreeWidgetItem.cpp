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
    int sortCol = treeWidget()->sortColumn();
    int myNumber = text(sortCol).toInt();
    int otherNumber = other.text(sortCol).toInt();
    return myNumber < otherNumber;
}