#include <iostream>
#include <Eigen/Core>
//#include <QGuiApplication>
#include <Qt3DCore>
#include <Qt3DRender>
#include <Qt3DExtras>


#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtWidgets/QHBoxLayout>
#include <MolecularGeometry3D.h>

#include "ElementInfo.h"
#include "ElementTypes.h"
#include "Sphere.h"
#include "Cylinder.h"
#include "DividedCylinder.h"
#include "Atom3D.h"
#include "Bond3D.h"
#include "Electron3D.h"
#include "Polyline.h"

#include "WaveFunctionParser.h"
#include "Atom.h"

#include "ArcLengthParametrizedBSpline.h"
#include "StringMethodCoordinatesPlotter.h"
#include "StringMethodValuePlotter.h"
#include "StringMethod.h"


int main(int argc, char *argv[]) {
  //Electron3D e1(root, QVector3D(12.0f, 1.0f, 7.0f), Spin::SpinType::Alpha);
  //Electron3D e2(root, QVector3D(9.0f, 1.2f, 7.0f), Spin::SpinType::Beta);

  //Qt3DRender::QObjectPicker(); // emits signals for you to handle
  //Qt3DRender::QPickingSettings* pickingSettings = new Qt3DRender::QPickingSettings();
  //pickingSettings->setPickMethod(Qt3DRender::QPickingSettings::PickMethod::BoundingVolumePicking);

  //QGuiApplication app(argc, argv);
  QApplication app(argc, argv);

  Eigen::VectorXd xA(2 * 3);
  Eigen::VectorXd xB(2 * 3);
  xA << 0.0, 0.0, +0.70014273, 0.0, 0.00, -0.70014273;
  xB << 0.0, 0.0, -0.70014273, 0.0, 0.00, -0.70014273;

  unsigned numberOfStates = 20;
  Eigen::MatrixXd initialCoordinates(2 * 3, numberOfStates);

  Eigen::VectorXd delta = xB - xA;
  for (int i = 0; i < numberOfStates; ++i) {
    double rel = double(i) / double(numberOfStates - 1);

    Eigen::VectorXd randVec(2*3);
    randVec.setZero();
    if ( i != 0  && i != (numberOfStates-1)) randVec.setRandom();
    randVec *= 0.25;

    initialCoordinates.col(i) = xA + ((delta+randVec) * rel);
  }

  StringMethod stringMethod(initialCoordinates);
  std::cout << stringMethod.getChain().coordinates() << std::endl;
  stringMethod.optimizeString();

  Qt3DCore::QEntity *root = new Qt3DCore::QEntity();

  // draw molecular geometry
  std::string filename = "t.wf";
  WaveFunctionParser waveFunctionParser(filename);
  waveFunctionParser.readNuclei();
  MolecularGeometry3D molecularGeometry3D (root, waveFunctionParser.getAtomCollection());


  BSplines::ArcLengthParametrizedBSpline bs = stringMethod.getArcLengthParametrizedBSpline();
  StringMethodCoordinatesPlotter bSplinePlotter(root,bs, 50, 0.005f);

  Qt3DExtras::Qt3DWindow *view = new Qt3DExtras::Qt3DWindow();
  view->defaultFrameGraph()->setClearColor(Qt::gray);

  QWidget *container = QWidget::createWindowContainer(view);

  Qt3DCore::QEntity *scene = root;

  // camera
  Qt3DRender::QCamera *camera = view->camera();
  camera->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
  camera->setPosition(QVector3D(0, 0, 2.0));
  camera->setViewCenter(QVector3D(0, 0, 0));
  //camera->setFarPlane(1000);
  //camera->setNearPlane(0.1);

  // manipulator
  //Qt3DExtras::QFirstPersonCameraController* manipulator = new Qt3DExtras::QFirstPersonCameraController(scene);
  Qt3DExtras::QOrbitCameraController *manipulator = new Qt3DExtras::QOrbitCameraController(scene);
  manipulator->setLinearSpeed(50.f);
  manipulator->setLookSpeed(180.f);
  manipulator->setCamera(camera);

  view->setRootEntity(scene);

  StringMethodValuesPlotter stringMethodValuesPlotter;
  QtCharts::QLineSeries *series = stringMethodValuesPlotter.getLineSeries(bs,200);

  QtCharts::QChart *chart = new QtCharts::QChart();
  chart->legend()->hide();
  chart->addSeries(series);
  chart->createDefaultAxes();
  chart->setTitle("Probability Density");

  QtCharts::QChartView *chartView = new QtCharts::QChartView(chart);
  chartView->setRenderHint(QPainter::Antialiasing);

  /*
  QMainWindow window;
  window.setCentralWidget(chartView);
  window.resize(400, 300);
  window.show();
  */

  QWidget *widget = new QWidget;
  QVBoxLayout *vLayout = new QVBoxLayout(widget);
  vLayout->addWidget(container,Qt::AlignTop);
  vLayout->addWidget(chartView,Qt::AlignBottom);

  widget->setWindowTitle(QStringLiteral("String Method Result"));

  widget->show();
  widget->resize(800, 600);

  return app.exec();

}
