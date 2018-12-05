//
// Created by heuer on 05.12.18.
//

#ifndef INPSIGHTS_MAXIMATREEWIDGETITEM_H
#define INPSIGHTS_MAXIMATREEWIDGETITEM_H

#include <QTreeWidget>

class MaximaTreeWidgetItem : public QTreeWidgetItem {
public:
    MaximaTreeWidgetItem(QTreeWidget *tree);

    MaximaTreeWidgetItem(QTreeWidgetItem *parent, const QStringList &strings, int type = Type);

    MaximaTreeWidgetItem(QTreeWidget * parent, const QStringList & strings);

    bool operator< (const QTreeWidgetItem &other) const override;
};

#endif //INPSIGHTS_MAXIMATREEWIDGETITEM_H
