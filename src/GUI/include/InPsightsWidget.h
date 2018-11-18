//
// Created by Michael Heuer on 18.11.18.
//

#ifndef INPSIGHTS_INPSIGHTSWIDGET_H
#define INPSIGHTS_INPSIGHTSWIDGET_H

#include <QWidget>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QSpinBox>
#include <MoleculeWidget.h>
#include <QTimer>
#include <QSplashScreen>

class InPsightsWidget : public QWidget{
    Q_OBJECT
public:

    explicit InPsightsWidget(QWidget *parent = nullptr)
    :
    QWidget(parent),
    layout_(new QHBoxLayout(this)),
    moleculeWidget_(new MoleculeWidget()),
    integerSpinBox_(new QSpinBox())
    {

        auto splashScreen = createSplashScreen();

        layout_->addWidget(moleculeWidget_, Qt::AlignLeft);
        layout_->addWidget(integerSpinBox_);

        integerSpinBox_->setRange(0, 0);
        integerSpinBox_->setSingleStep(1);
        integerSpinBox_->setValue(0);

        setLayout(layout_);
        resize(800,800);

        QTimer::singleShot(2000, splashScreen, SLOT(close()));
        QTimer::singleShot(2000, this, SLOT(show()));
    }

    QHBoxLayout *layout_;
    MoleculeWidget *moleculeWidget_;
    QSpinBox *integerSpinBox_;


private:
    QSplashScreen* createSplashScreen(){
        auto splashScreen = new QSplashScreen();

        QImage file(":inPsights.png");
        QPixmap qPixmap;
        qPixmap.convertFromImage(file);
        splashScreen->setPixmap(qPixmap);
        splashScreen->show();
        splashScreen->showMessage("Version 1.0.0", Qt::AlignRight, Qt::lightGray);

        return splashScreen;
    }
};

#endif //INPSIGHTS_INPSIGHTSWIDGET_H
