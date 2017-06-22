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
#include "StationaryPointFinder.h"

#include "solver/newtonraphsonsolver.h"
#include "ElectronicWaveFunctionProblem.h"
#include <Eigen/Eigenvalues>


int main(int argc, char *argv[]) {
  //Electron3D e1(root, QVector3D(12.0f, 1.0f, 7.0f), Spin::SpinType::Alpha);
  //Electron3D e2(root, QVector3D(9.0f, 1.2f, 7.0f), Spin::SpinType::Beta);

  //Qt3DRender::QObjectPicker(); // emits signals for you to handle
  //Qt3DRender::QPickingSettings* pickingSettings = new Qt3DRender::QPickingSettings();
  //pickingSettings->setPickMethod(Qt3DRender::QPickingSettings::PickMethod::BoundingVolumePicking);

  //QGuiApplication app(argc, argv);
  QApplication app(argc, argv);


  // H2
  //Eigen::VectorXd xA(2 * 3);
  //Eigen::VectorXd xB(2 * 3);
  //xA << 0.0, 0.0, +0.70014273, 0.0, 0.00, -0.70014273;
  //xB << 0.0, 0.0, -0.70014273, 0.0, 0.00, -0.70014273;

  Eigen::VectorXd xA(18*3);
  Eigen::VectorXd xB(18*3);
  /*CTR*/
//  xA <<  \
// 0.000000, 0.000000, 1.443089,\
// 0.000000, 0.000000,-1.443089,\
//-1.663110, 0.960208,-2.193991,\
// 0.000000, 1.920396, 2.193991,\
//-1.663110,-0.960208, 2.193991,\
// 0.000637,-0.881335,-1.752072,\
// 0.762826, 0.441312,-1.754787,\
// 0.762219,-0.441827, 1.754450,\
// 0.000696, 0.003032, 0.695882,/*9*/\
// 0.000000, 0.000000, 1.443089,\
// 0.000000, 0.000000,-1.443089,\
// 1.663110,-0.960208, 2.193991,\
// 0.000000,-1.920396,-2.193991,\
// 1.663110, 0.960208,-2.193991,\
//-0.000261, 0.881214, 1.752389,\
//-0.762453,-0.441670, 1.755289,\
//-0.001111,-0.002637,-0.695660,\
//-0.762504, 0.441871,-1.753631/*18*/;

  /*max*/
  xA << \
 0.000000, 0.000000, 1.443089,\
 0.000000, 0.000000,-1.443089,\
-1.663110, 0.960208,-2.193991,\
 0.000000, 1.920396, 2.193991,\
-1.663110,-0.960208, 2.193991,\
-0.023346,-0.882838,-1.743612,\
 0.752879, 0.461638,-1.743608,\
 0.756924,-0.437014, 1.772883,\
-0.044372, 0.025617, 0.695896,\
 0.000000, 0.000000, 1.443089,\
 0.000000, 0.000000,-1.443089,\
 1.663110,-0.960208, 2.193991,\
 0.000000,-1.920396,-2.193991,\
 1.663110, 0.960208,-2.193991,\
 0.023348, 0.882821, 1.743604,\
-0.752862,-0.461629, 1.743600,\
 0.044374,-0.025617,-0.695900,\
-0.756913, 0.437007,-1.772876;

  Eigen::Vector3d el9= xA.segment((9-1)*3,3);
  Eigen::Vector3d el17= xA.segment((17-1)*3,3);

  //Eigen::Vector3d el6= xA.segment((6-1)*3,3);
  //Eigen::Vector3d el7= xA.segment((7-1)*3,3);
  //Eigen::Vector3d el18= xA.segment((18-1)*3,3);


  xB = xA;
  xB.segment((9-1)*3,3) = el17;
  xB.segment((17-1)*3,3) = el9;
  //xB.segment((6-1)*3,3) = el7;
  //xB.segment((7-1)*3,3) = el18;
  //xB.segment((18-1)*3,3) = el6;

  ElectronicWaveFunction::getInstance().evaluate(xA);
  std::cout << "phi " << ElectronicWaveFunction::getInstance().getDeterminantProbabilityAmplitude() << std::endl;
  std::cout << "U " << ElectronicWaveFunction::getInstance().getJastrowFactor() << std::endl;
  std::cout << "phi*exp(U) " << ElectronicWaveFunction::getInstance().getProbabilityAmplitude() << std::endl;
  std::cout << "eloc " << ElectronicWaveFunction::getInstance().getLocalEnergy() << std::endl;
  std::cout << "drift " << ElectronicWaveFunction::getInstance().getElectronDriftCollection().transpose() << std::endl;

  ElectronicWaveFunction::getInstance().evaluate(xB);
  std::cout << ElectronicWaveFunction::getInstance().getDeterminantProbabilityAmplitude() << std::endl;


  unsigned numberOfStates = 10;
  Eigen::MatrixXd initialCoordinates(18 * 3, numberOfStates);
  Eigen::VectorXd delta = xB - xA;
  for (int i = 0; i < numberOfStates; ++i) {
    double rel = double(i) / double(numberOfStates - 1);
    initialCoordinates.col(i) = xA + (delta * rel);
  }

  //Eigen::VectorXd bend(numberOfStates);
  //bend << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  //x-bend
  //initialCoordinates.row((9-1)*3+0) += bend; //+x bend
  //initialCoordinates.row((17-1)*3+0) -= bend;//-x bend

  //y-bend
  //initialCoordinates.row((9-1)*3+1) += bend; //+y bend
  //initialCoordinates.row((17-1)*3+1) -= bend;//-y bend

  StringMethod stringMethod(initialCoordinates);
  stringMethod.optimizeString();

  BSplines::ArcLengthParametrizedBSpline arcLengthParametrizedBSpline = stringMethod.getArcLengthParametrizedBSpline();
  BSplines::BSpline bspline(arcLengthParametrizedBSpline.getKnotVector(0),
                       arcLengthParametrizedBSpline.getControlPointMatrix(0),
                       arcLengthParametrizedBSpline.getDegree(),2);

  BSplines::StationaryPointFinder stationaryPointFinder(bspline);
  std::vector<double> result = stationaryPointFinder.getMaxima(0);

  Eigen::VectorXd tsGuessGeom = bspline.evaluate(0).tail(18*3);
  std::cout << "u=" << result[0] << "\n" << tsGuessGeom.transpose() << std::endl;


  ElectronicWaveFunctionProblem f;
  Eigen::MatrixXd hess(18*3,18*3);
  f.hessian(tsGuessGeom,hess);
  //std::cout << hess << std::endl;
  Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(hess,false);
  auto eigenvalues = eigenSolver.eigenvalues();
  std::cout << eigenvalues << std::endl;


    //Eigen::VectorXd tsGuessGeom = arcLengthParametrizedBSpline.evaluate(*it).tail(18*3);

    cppoptlib::NewtonRaphsonSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
    crit.iterations = 100;
    crit.gradNorm = 1e-3;
    solver.setStopCriteria(crit);
    solver.minimize(f,tsGuessGeom);

  f.hessian(tsGuessGeom,hess);
  //std::cout << hess << std::endl;
  eigenSolver = Eigen::EigenSolver<Eigen::MatrixXd>(hess,false);
  eigenvalues = eigenSolver.eigenvalues();
  std::cout << eigenvalues << std::endl;





  Qt3DCore::QEntity *root = new Qt3DCore::QEntity();

  // draw molecular geometry
  std::string filename = "t.wf";
  WaveFunctionParser waveFunctionParser(filename);
  waveFunctionParser.readNuclei();
  MolecularGeometry3D molecularGeometry3D (root, waveFunctionParser.getAtomCollection());

  for (int j = 0; j < xA.rows()/3 ; ++j) {
    Eigen::Vector3d vec3d = xA.segment(j*3,3);
    QVector3D qVector3D (vec3d(0),vec3d(1),vec3d(2));

    if ( j < (xA.rows()/3)/2 )
      Electron3D(root,qVector3D,Spin::Alpha);
    else
      Electron3D(root,qVector3D,Spin::Beta);
  }

  BSplines::ArcLengthParametrizedBSpline bs = stringMethod.getArcLengthParametrizedBSpline();
  StringMethodCoordinatesPlotter bSplinePlotter(root,bs, 50, 0.005f);

  Qt3DExtras::Qt3DWindow *view = new Qt3DExtras::Qt3DWindow();
  view->defaultFrameGraph()->setClearColor(Qt::gray);

  QWidget *container = QWidget::createWindowContainer(view);

  Qt3DCore::QEntity *scene = root;

  // camera
  Qt3DRender::QCamera *camera = view->camera();
  camera->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
  camera->setPosition(QVector3D(6.0, 0.0, 0.0));
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


  QWidget *widget = new QWidget;
  QVBoxLayout *vLayout = new QVBoxLayout(widget);
  vLayout->addWidget(container,2);
  vLayout->addWidget(chartView,1);

  widget->setWindowTitle(QStringLiteral("String Method Result"));

  widget->show();
  widget->resize(800, 800);

  return app.exec();

}
