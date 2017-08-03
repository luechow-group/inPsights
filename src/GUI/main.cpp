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
#include <QtCharts/QValueAxis>
#include <QtWidgets/QHBoxLayout>
#include <Qt3DExtras>

#include "MolecularGeometry3D.h"

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
#include <QtCharts/QScatterSeries>


int main(int argc, char *argv[]) {
  //Qt3DRender::QObjectPicker(); // emits signals for you to handle
  //Qt3DRender::QPickingSettings* pickingSettings = new Qt3DRender::QPickingSettings();
  //pickingSettings->setPickMethod(Qt3DRender::QPickingSettings::PickMethod::BoundingVolumePicking);

  //QGuiApplication app(argc, argv);
  QApplication app(argc, argv);

  ElectronicWaveFunction::getInstance("Ethane-em-5.wf");
  ElectronicWaveFunction::getInstance().setFrozenElectrons({1,2,3,4,5, 10,11,12,13,14});

  Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
  xA << \
 0.000000, 0.000000, 1.443184,\
 0.000000, 0.000000,-1.443184,\
-1.662146,-0.959641, 2.192989,\
-1.662146, 0.959641,-2.192989,\
 0.000000,-1.919300,-2.192989,\
-0.024099, 0.773535, 1.718336,\
 0.657845,-0.407636, 1.718335,\
 0.658424, 0.380140,-1.752435,\
-0.034848,-0.020119,-0.660100,\
 0.000000, 0.000000, 1.443184,\
 0.000000, 0.000000,-1.443184,\
 1.662146,-0.959641, 2.192989,\
 1.662146, 0.959641,-2.192989,\
 0.000000, 1.919300, 2.192989,\
-0.657846, 0.407635,-1.718336,\
 0.024100,-0.773537,-1.718336,\
-0.658423,-0.380140, 1.752435,\
 0.034847, 0.020120, 0.660096;

  Eigen::Vector3d el8= xA.segment((8-1)*3,3);
  Eigen::Vector3d el9= xA.segment((9-1)*3,3);
  Eigen::Vector3d el17= xA.segment((17-1)*3,3);
  Eigen::Vector3d el18= xA.segment((18-1)*3,3);

  Eigen::Vector3d el15= xA.segment((15-1)*3,3);
  Eigen::Vector3d el16= xA.segment((16-1)*3,3);

  Eigen::VectorXd xB(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
  xB = xA;
  //xB.segment((9-1)*3,3) = el18;
  //xB.segment((18-1)*3,3) = el9;
  xB.segment((9-1)*3,3) = el8;
  xB.segment((8-1)*3,3) = el9;
  //xB.segment((17-1)*3,3) = el18;
  //xB.segment((18-1)*3,3) = el17;

/*
  ElectronicWaveFunction::getInstance("Ethylene-em-5.wf");
  ElectronicWaveFunction::getInstance().setFrozenElectrons({1,2,3,4, 9,10,11,12});

  Eigen::VectorXd xA(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
  xA <<\
      0.00000000, 0.00000000, 1.25132000,\
 0.00000000, 0.00000000,-1.25132000,\
 0.00000000, 1.73896400,-2.32682000,\
 0.00000000, 1.73896400, 2.32682000,\
 0.06517605,-0.70088655,-1.73499063,\
-0.06517605,-0.70088655, 1.73499063,\
-0.56040160, 0.00175003,-0.79662456,\
 0.56040160, 0.00175003, 0.79662456,\
 0.00000000, 0.00000000, 1.25132000,\
 0.00000000, 0.00000000,-1.25132000,\
 0.00000000,-1.73896400, 2.32682000,\
 0.00000000,-1.73896400,-2.32682000,\
-0.06517605, 0.70088655,-1.73499063,\
 0.06517605, 0.70088655, 1.73499063,\
 0.56040160,-0.00175003,-0.79662456,\
-0.56040160,-0.00175003, 0.79662456;


  Eigen::VectorXd xB(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
  xB <<\
 0.00000000, 0.00000000, 1.25132000,\
 0.00000000, 0.00000000,-1.25132000,\
 0.00000000, 1.73896400,-2.32682000,\
 0.00000000, 1.73896400, 2.32682000,\
 0.00001705,-0.68411983,-1.75563538,\
 0.00000116,-0.72989937, 1.71567587,\
-0.55865730, 0.01191069,-0.79672667,\
 0.55863148, 0.01192229,-0.79668313,\
 0.00000000, 0.00000000, 1.25132000,\
 0.00000000, 0.00000000,-1.25132000,\
 0.00000000,-1.73896400, 2.32682000,\
 0.00000000,-1.73896400,-2.32682000,\
 0.00000116, 0.72989937,-1.71567587,\
 0.00001705, 0.68411983, 1.75563538,\
 0.55863148,-0.01192229, 0.79668313,\
-0.55865730,-0.01191069, 0.79672667;*/


/*
  // Direct rectangle exchange
  Eigen::Vector3d el7= xA.segment((7-1)*3,3);
  Eigen::Vector3d el8= xA.segment((8-1)*3,3);
  Eigen::Vector3d el15= xA.segment((15-1)*3,3);
  Eigen::Vector3d el16= xA.segment((16-1)*3,3);
  Eigen::VectorXd xB(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);

  xB = xA;
  xB.segment((7-1)*3,3)  = el8;
  xB.segment((8-1)*3,3)  = el7;
  xB.segment((15-1)*3,3) = el16;
  xB.segment((16-1)*3,3) = el15;*/


  ElectronicWaveFunction::getInstance().evaluate(xA);
  std::cout << "phi " << ElectronicWaveFunction::getInstance().getDeterminantProbabilityAmplitude() << std::endl;
  std::cout << "U " << ElectronicWaveFunction::getInstance().getJastrowFactor() << std::endl;
  std::cout << "phi*exp(U) " << ElectronicWaveFunction::getInstance().getProbabilityAmplitude() << std::endl;
  std::cout << "eloc " << ElectronicWaveFunction::getInstance().getLocalEnergy() << std::endl;
  std::cout << "drift " << ElectronicWaveFunction::getInstance().getElectronDriftCollection().transpose() << std::endl;

  ElectronicWaveFunction::getInstance().evaluate(xB);
  std::cout << ElectronicWaveFunction::getInstance().getDeterminantProbabilityAmplitude() << std::endl;


  unsigned numberOfStates = 11;

  Eigen::MatrixXd initialCoordinates(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3, numberOfStates);
  Eigen::VectorXd delta = xB - xA;
  for (int i = 0; i < numberOfStates; ++i) {
    double rel = double(i) / double(numberOfStates - 1);
    initialCoordinates.col(i) = xA + (delta * rel);
  }

  /*Eigen::VectorXd bend1(numberOfStates);
  bend1  << 0,0.1,0.3,0.4,0.5,0.6,0.8,0.6,0.3,0.1,0;
  initialCoordinates.row((15-1)*3+1) -= bend1; //+y bend

  Eigen::VectorXd bend2(numberOfStates);
  bend2  << 0,0,0.2,0.4,0.3,0.2,0.1,0,0,0,0;
  initialCoordinates.row((8-1)*3+1) += 0.5*bend2; //+y bend*/

  // BEND FUNCTION
  double shift = 0.4;
  Eigen::VectorXd bend(numberOfStates);
  for (int k = 0; k < numberOfStates/2+1; ++k) {
    bend(k) = std::sqrt(shift/(numberOfStates/2)*k);
  }
  //bend(numberOfStates/2+1) = bend(numberOfStates/2);

  for (int k = numberOfStates-1; k >= numberOfStates/2; --k) {
    bend(k) = bend(numberOfStates-k-1);
  }
  std::cout << bend.transpose() << std::endl;

  //y-bend
  initialCoordinates.row((8-1)*3+1) -= 0.5*bend;//-y bend
  initialCoordinates.row((9-1)*3+1) += 0.5*bend; //+y bend
  //initialCoordinates.row((17-1)*3+1) -= 0.5*bend;//-y bend
  //initialCoordinates.row((18-1)*3+1) += 0.5*bend; //+y bend



  StringMethod stringMethod(initialCoordinates);
  stringMethod.optimizeString();

  BSplines::ArcLengthParametrizedBSpline arcLengthParametrizedBSpline = stringMethod.getArcLengthParametrizedBSpline();
  BSplines::BSpline bspline(arcLengthParametrizedBSpline.getKnotVector(0),
                       arcLengthParametrizedBSpline.getControlPointMatrix(0),
                       arcLengthParametrizedBSpline.getDegree(),2);


  BSplines::StationaryPointFinder stationaryPointFinder(bspline);
  std::vector<double> minima = stationaryPointFinder.getMinima(0);
  std::vector<double> maxima = stationaryPointFinder.getMaxima(0);

  for (int l = 0; l < minima.size(); ++l) {
    std::cout << "l=" << l << ", u="<< minima[l]<< std::endl;
    std::cout << bspline.evaluate(minima[l]).tail(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3).transpose() << std::endl;
  }

  std::cout << bspline.evaluate(0).tail(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3).transpose() << std::endl;
  std::cout << bspline.evaluate(1).tail(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3).transpose() << std::endl;

  Eigen::VectorXd tsGuessGeom = bspline.evaluate(maxima[0]).tail(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
  std::cout << "u=" << maxima[0] << "\n" << tsGuessGeom.transpose() << std::endl;
  std::cout << "0_to_ts: " << bspline.evaluate(maxima[0])(0)-bspline.evaluate(0)(0) << std::endl;
  std::cout << "ts_to_1: " << bspline.evaluate(maxima[0])(0)-bspline.evaluate(1)(0) << std::endl;

  /*
  ElectronicWaveFunctionProblem f("");
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
    crit.iterations = 15;
    crit.gradNorm = 1e-5;
    solver.setStopCriteria(crit);
    solver.minimize(f,tsGuessGeom);

  f.hessian(tsGuessGeom,hess);
  //std::cout << hess << std::endl;
  eigenSolver = Eigen::EigenSolver<Eigen::MatrixXd>(hess,false);
  eigenvalues = eigenSolver.eigenvalues();
  std::cout << eigenvalues << std::endl;*/


  Qt3DCore::QEntity *root = new Qt3DCore::QEntity();

  //Draw tsguess
  for (int j = 0; j < ElectronicWaveFunction::getInstance().getNumberOfElectrons() ; ++j) {
    Eigen::Vector3d vec3d = tsGuessGeom.segment(j * 3, 3);
    QVector3D qVector3D(vec3d(0), vec3d(1), vec3d(2));

    Electron3D* e = new Electron3D(root, qVector3D, Spin::None);
    e->setRadius(0.025f);
  }

  // draw molecular geometry
  std::cout << "wf:" << ElectronicWaveFunction::getInstance().getFileName() << std::endl;
  WaveFunctionParser waveFunctionParser(ElectronicWaveFunction::getInstance().getFileName());
  waveFunctionParser.readNuclei();
  MolecularGeometry3D molecularGeometry3D (root, waveFunctionParser.getAtomCollection());

  // draw electrons
  for (int j = 0; j < ElectronicWaveFunction::getInstance().getNumberOfElectrons() ; ++j) {
    Eigen::Vector3d vec3d = xA.segment(j*3,3);
    QVector3D qVector3D (vec3d(0),vec3d(1),vec3d(2));

    Electron3D *e;
    if ( j < (xA.rows()/3)/2 ){
      e = new Electron3D(root,qVector3D,Spin::Alpha);
    }
    else {
      e = new Electron3D(root, qVector3D, Spin::Beta);
    }

    /*
    auto *textMaterial = new Qt3DExtras::QPhongMaterial(e);
    { // text
        auto *text = new Qt3DCore::QEntity(e);
        auto *textMesh = new Qt3DExtras::QExtrudedTextMesh();

        auto *textTransform = new Qt3DCore::QTransform();
        QFont font(QString("Arial"), 12
                , 0, false);
      QVector3D shift;

      if(e->getSpinType() == Spin::Alpha) shift = QVector3D(0.0f,0.07f,0.07f);
      else shift = QVector3D(-0.0f,-0.07f,-0.07f);

        textTransform->setTranslation(qVector3D+shift);
        textTransform->setRotationY(90);
        textTransform->setScale(.05f);
        textMesh->setDepth(.1f);
        textMesh->setFont(font);
        textMesh->setText(QString::fromStdString(std::to_string(j+1)));
        textMaterial->setDiffuse(e->getColor());

        text->addComponent(textMaterial);
        text->addComponent(textMesh);
        text->addComponent(textTransform);
    }
     */
  }

  BSplines::ArcLengthParametrizedBSpline bs = stringMethod.getArcLengthParametrizedBSpline();
  StringMethodCoordinatesPlotter bSplinePlotter(root,bs, 100, 0.01f);

  Qt3DExtras::Qt3DWindow *view = new Qt3DExtras::Qt3DWindow();
  view->defaultFrameGraph()->setClearColor(Qt::white);

  QWidget *moleculeView = QWidget::createWindowContainer(view);

  Qt3DCore::QEntity *scene = root;
  // camera
  Qt3DRender::QCamera *camera = view->camera();
  camera->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
  //camera->setPosition(QVector3D(4.25, 4.25, 0.0)); // ethane
  camera->setPosition(QVector3D(2.5, -5.00, 0.0)); // ethane
  camera->setViewCenter(QVector3D(0, 0, 0));
  //camera->setFarPlane(1000);
  //camera->setNearPlane(0.1);
  //camera->tilt(-5.0f);

  // manipulator
  //Qt3DExtras::QFirstPersonCameraController* manipulator = new Qt3DExtras::QFirstPersonCameraController(scene);
  Qt3DExtras::QOrbitCameraController *manipulator = new Qt3DExtras::QOrbitCameraController(scene);
  manipulator->setLinearSpeed(50.f);
  manipulator->setLookSpeed(180.f);
  manipulator->setCamera(camera);


  view->setRootEntity(scene);

  StringMethodValuesPlotter stringMethodValuesPlotter;
  QtCharts::QLineSeries *lineSeries = stringMethodValuesPlotter.getLineSeries(bs,200);
  lineSeries->setColor(Qt::black);
  lineSeries->setName("string");

  QtCharts::QScatterSeries *scatterSeries = new QtCharts::QScatterSeries();
  scatterSeries->setColor(Qt::gray);
  scatterSeries->setName("states");
  scatterSeries->setMarkerSize(10.0);
  Eigen::VectorXd values = stringMethod.getChain().values();
  for (int i = 0; i < values.size(); ++i) {

    double u = double(i)/double(values.size());

    scatterSeries->append(u,values(i));
  }



  QtCharts::QChartView *chartView = new QtCharts::QChartView();
  chartView->setRenderHint(QPainter::Antialiasing);

  //chartView->chart()->legend()->hide();
  chartView->chart()->addSeries(lineSeries);
  chartView->chart()->addSeries(scatterSeries);
  //chartView->chart()->createDefaultAxes();
  //chartView->chart()->setTitle("Probability Density");

  QtCharts::QValueAxis *axisX = new QtCharts::QValueAxis;
  axisX->setTitleText("rel. (weighted) arc length <i>u</i>");
  //axisX->setRange(0, 1);
  //axisX->setTickCount(10);
  chartView->chart()->setAxisX(axisX,lineSeries);

  QtCharts::QValueAxis *axisY = new QtCharts::QValueAxis;
  axisY->setTitleText("-ln(|Ψ|²)");
  chartView->chart()->setAxisY(axisY,lineSeries);


  QWidget *widget = new QWidget;
  QVBoxLayout *vLayout = new QVBoxLayout(widget);
  vLayout->addWidget(moleculeView,2);
  vLayout->addWidget(chartView,1);

  widget->setWindowTitle(QStringLiteral("String Method Result"));

  widget->show();
  widget->resize(1920, 1080);

  return app.exec();

}
