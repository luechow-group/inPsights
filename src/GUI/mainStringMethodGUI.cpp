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

#include "MoleculeWidget.h"
#include "AtomCollection3D.h"
#include "ElectronCollection3D.h"

#include "ElementInfo.h"
#include "ElementType.h"
#include "Sphere.h"
#include "Bond3D.h"
#include "Electron3D.h"
#include "Polyline.h"

#include "WfFileImporter.h"
#include "Atom.h"

#include "ArcLengthParametrizedBSpline.h"
#include "StringMethodCoordinatesPlotter.h"
#include "StringMethodValuePlotter.h"
#include "StringMethod.h"
#include "StationaryPointFinder.h"

#include "ElectronicWaveFunctionProblem.h"
#include <Eigen/Eigenvalues>
#include <QtCharts/QScatterSeries>

#include "solver/bfgsnssolver.h"
#include "CollectionParser.h"

int main(int argc, char *argv[]) {

    QApplication app(argc, argv);

    //ElectronicWaveFunction::getInstance("H2ic666.wf");//"Ethane-em-5.wf");
    ElectronicWaveFunctionProblem f("H2ic666.wf");
    std::cout << f.getNucleiPositions() << std::endl;
    Eigen::VectorXd x0(2*3);
    x0 <<
       0,0,-0.700144,\
       0,0,+0.700144;



    /*
    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();
    crit.iterations = 1000;
    crit.gradNorm = 1e-8;
    cppoptlib::BfgsnsSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    solver.setStopCriteria(crit);
    Eigen::VectorXd x = x0;
    solver.minimize(f, x);
    std::cout << "max: " << x << std::endl;*/



    //collectionParser.writeJSON(collectionParser.electronCollectionToJson(ec),"Ethylene-glob-max.json");
    CollectionParser collectionParser;
    //auto ecA = collectionParser.electronCollectionFromJson("Ethane-glob-max.json");
    auto ecA = ElectronCollection(x0);
    auto ecB = ecA;
    ecB.permute(0, 1);
    //ecB.permute(8, 17);//Ethane

    unsigned numberOfStates = 7;

    auto xA = ecA.positionsAsEigenVector();
    auto xB = ecB.positionsAsEigenVector();

    Eigen::MatrixXd initialCoordinates(ElectronicWaveFunction::getInstance().getNumberOfElectrons() * 3,
                                       numberOfStates);
    Eigen::VectorXd delta = xB - xA;
    for (int i = 0; i < numberOfStates; ++i) {
        double rel = double(i) / double(numberOfStates - 1);
        initialCoordinates.col(i) = xA + (delta * rel);
    }


    Eigen::VectorXd bend(numberOfStates);
    bend << 0, 0.0, 0.05, 0.1, \
            //0.15, 0.1,
            0.05, 0.0, 0;
    //bend.reverse();
    std::cout << bend << std::endl;
    initialCoordinates.row(0*3 + 0) += 0.4 * bend;//x bend
    initialCoordinates.row(1*3 + 0) -= 0.4 * bend;
    initialCoordinates.row(0*3 + 1) += 0.4 * bend;//y bend
    initialCoordinates.row(1*3 + 1) -= 0.4 * bend;


    /*Ethanae
Eigen::VectorXd bend(numberOfStates);
bend << 0, 0.0, 0.05, 0.10, 0.3, 0.5, 0.4, 0.1, 0;
//bend.reverse();
std::cout << bend << std::endl;

//Ethane Bend C-C bond exchange
initialCoordinates.row((9 - 1) * 3 + 0) -= 0.1 * bend.reverse();//x bend
initialCoordinates.row((18 - 1) * 3 + 0) += 0.1 * bend; //-x bend
initialCoordinates.row((9 - 1) * 3 + 1) += 0.05 * bend;//y bend
initialCoordinates.row((18 - 1) * 3 + 2) -= 0.05 * bend;//z bend
*/

    StringMethod stringMethod(initialCoordinates);
    stringMethod.optimizeString();

    BSplines::ArcLengthParametrizedBSpline arcLengthParametrizedBSpline = stringMethod.getArcLengthParametrizedBSpline();
    BSplines::BSpline bspline(arcLengthParametrizedBSpline.getKnotVector(0),
                              arcLengthParametrizedBSpline.getControlPointMatrix(0),
                              arcLengthParametrizedBSpline.getDegree(), 2);


    BSplines::StationaryPointFinder stationaryPointFinder(bspline);
    std::vector<double> minima = stationaryPointFinder.getMinima(0);
    std::vector<double> maxima = stationaryPointFinder.getMaxima(0);

    for (int l = 0; l < minima.size(); ++l) {
        std::cout << "l=" << l << ", u=" << minima[l] << std::endl;
        std::cout << bspline.evaluate(minima[l]).tail(
                ElectronicWaveFunction::getInstance().getNumberOfElectrons() * 3).transpose() << std::endl;
    }
    Eigen::VectorXd xAopt(
            bspline.evaluate(0).tail(ElectronicWaveFunction::getInstance().getNumberOfElectrons() * 3).transpose());
    Eigen::VectorXd xBopt(
            bspline.evaluate(1).tail(ElectronicWaveFunction::getInstance().getNumberOfElectrons() * 3).transpose());

    std::cout << "u=0\n" << xAopt.transpose() << std::endl;
    std::cout << "u=1\n" << xBopt.transpose() << std::endl;

    Eigen::VectorXd tsGuessGeom = bspline.evaluate(maxima[0]).tail(
            ElectronicWaveFunction::getInstance().getNumberOfElectrons() * 3);
    std::cout << "u=" << maxima[0] << "\n" << tsGuessGeom.transpose() << std::endl;
    std::cout << "0_to_ts: " << bspline.evaluate(maxima[0])(0) - bspline.evaluate(0)(0) << std::endl;
    std::cout << "ts_to_1: " << bspline.evaluate(maxima[0])(0) - bspline.evaluate(1)(0) << std::endl;


    /* Eigen::MatrixXd hess(ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3,ElectronicWaveFunction::getInstance().getNumberOfElectrons()*3);
     f.hessian(xA,hess);
     //std::cout << hess << std::endl;
     Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(hess,false);
     auto eigenvalues = eigenSolver.eigenvalues();
     std::cout << eigenvalues << std::endl;

       //Eigen::VectorXd tsGuessGeom = arcLengthParametrizedBSpline.evaluate(*it).tail(18*3);
       cppoptlib::NewtonRaphsonSolver<ElectronicWaveFunctionProblem> solver;
       solver.setDebug(cppoptlib::DebugLevel::High);
       cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
       crit.iterations = 20;
       crit.gradNorm = 1e-7;
       solver.setStopCriteria(crit);
       solver.minimize(f,xA);

      f.hessian(xA,hess);
      //std::cout << hess << std::endl;
      eigenSolver = Eigen::EigenSolver<Eigen::MatrixXd>(hess,false);
      eigenvalues = eigenSolver.eigenvalues();
      std::cout << eigenvalues << std::endl;*/

    Eigen::VectorXd grad(ElectronicWaveFunction::getInstance().getNumberOfElectrons() * 3);
    f.gradient(xAopt, grad);
    std::cout << "u=0\n" << grad.transpose() << std::endl;
    f.gradient(tsGuessGeom, grad);
    std::cout << "u=" << maxima[0] << "\n" << grad.transpose() << std::endl;
    f.gradient(xBopt, grad);
    std::cout << "u=1\n" << grad.transpose() << std::endl;


    //visualization
    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    // draw molecular geometry
    std::cout << "wf:" << ElectronicWaveFunction::getInstance().getFileName() << std::endl;
    WfFileImporter waveFunctionParser(ElectronicWaveFunction::getInstance().getFileName());

    AtomCollection3D molecularGeometry3D(root, waveFunctionParser.getAtomCollection());
    ElectronCollection3D(root, ElectronicWaveFunction::getInstance().getElectronPositionCollection(), true);


    //Draw tsguess
    /*for (int j = 0; j < ElectronicWaveFunction::getInstance().getNumberOfElectrons() ; ++j) {
      Eigen::Vector3d vec3d = tsGuessGeom.segment(j * 3, 3);
      QVector3D qVector3D(vec3d(0), vec3d(1), vec3d(2));

      Electron3D* e = new Electron3D(root, qVector3D, Spin::SpinType::none);
      e->setRadius(0.025f);
      e->setAlpha(1.0);
    }*/

    BSplines::ArcLengthParametrizedBSpline bs = stringMethod.getArcLengthParametrizedBSpline();
    StringMethodCoordinatesPlotter bSplinePlotter(root, bs, 100, 0.01f);

    Qt3DExtras::Qt3DWindow *view = new Qt3DExtras::Qt3DWindow();
    view->defaultFrameGraph()->setClearColor(Qt::white);

    QWidget *moleculeView = QWidget::createWindowContainer(view);

    // camera
    Qt3DRender::QCamera *camera = view->camera();
    camera->lens()->setPerspectiveProjection(45.0f, 16.0f / 9.0f, 0.1f, 100.0f);
    camera->setPosition(QVector3D(2.5, -5.00, 0.0)); // ethane
    camera->setViewCenter(QVector3D(0, 0, 0));

    // manipulator
    Qt3DExtras::QOrbitCameraController *manipulator = new Qt3DExtras::QOrbitCameraController(root);
    manipulator->setLinearSpeed(50.f);
    manipulator->setLookSpeed(180.f);
    manipulator->setCamera(camera);

    view->setRootEntity(root);

    StringMethodValuesPlotter stringMethodValuesPlotter;
    QtCharts::QLineSeries *lineSeries = stringMethodValuesPlotter.getLineSeries(bs, 200);
    lineSeries->setColor(Qt::black);
    lineSeries->setName("string");

    QtCharts::QScatterSeries *scatterSeries = new QtCharts::QScatterSeries();
    scatterSeries->setColor(Qt::gray);
    scatterSeries->setName("states");
    scatterSeries->setMarkerSize(10.0);
    Eigen::VectorXd values = stringMethod.getChain().values();
    for (int i = 0; i < values.size(); ++i) {

        double u = double(i) / double(values.size());

        scatterSeries->append(u, values(i));
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
    chartView->chart()->setAxisX(axisX, lineSeries);

    QtCharts::QValueAxis *axisY = new QtCharts::QValueAxis;
    axisY->setTitleText("-ln(|Ψ|²)");
    chartView->chart()->setAxisY(axisY, lineSeries);


    QWidget *widget = new QWidget;
    QVBoxLayout *vLayout = new QVBoxLayout(widget);
    vLayout->addWidget(moleculeView, 2);
    vLayout->addWidget(chartView, 1);

    widget->setWindowTitle(QStringLiteral("Amolqc++"));

    widget->show();
    widget->resize(800, 600);

    return app.exec();

}
