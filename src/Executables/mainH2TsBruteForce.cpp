//
// Created by Michael Heuer on 01.03.18.
//

#include <QApplication>
#include <solver/newtonraphsonsolver.h>

#include "Serialization.h"
#include "ElectronicWaveFunctionProblem.h"
#include "solver/bfgsnssolver.h"
#include "solver/bfgssolver.h"
#include "solver/timeintegrationsolver.h"
#include "solver/gradientdescentumrigarlimitedsteplength.h"
#include "solver/gradientdescentsolver.h"
#include "solver/gradientdescentsimplesolver.h"
#include "solver/timeintegrationumrigarsolver.h"
#include "solver/bfgsumrigarsolver.h"
#include "solver/bfgsnssolver.h"

#include "AtomCollection3D.h"
#include "ElectronCollection3D.h"
#include "ParticleCollectionPath3D.h"
#include "MoleculeWidget.h"


int main(int argc, char *argv[]) {
    std::string wavefunctionFilename = "H2sm444.wf";
    bool showGui = true;

    ElectronicWaveFunctionProblem electronicWaveFunctionProblem(wavefunctionFilename);


    cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::nonsmoothDefaults();


    Eigen::VectorXd bothRight(6);
    bothRight<< 0,0,0.700144,0,0,0.700144;
    Eigen::VectorXd bothLeft(6);
    bothLeft << 0,0,-0.700144,0,0-0.700144;
    Eigen::VectorXd secondOrderTS1(6);
    secondOrderTS1 << 0, 0.19589114963364730,0,0, -0.1958911496336473,0;
    Eigen::VectorXd secondOrderTS2(6);
    secondOrderTS2 << 0,-0.19589114963364730,0,0, 0.1958911496336473,0;

    int mz = 4;
    int my = 4;
    double maxz = 0.35;//0.9*0.700144;
    double maxy = 0.3;

    Eigen::VectorXd vec(6);
    Eigen::VectorXd grad(6);
    double x = 0;

    cppoptlib::NewtonRaphsonSolver<ElectronicWaveFunctionProblem> solver;
    solver.setDebug(cppoptlib::DebugLevel::High);
    crit.gradNorm = 1e-6;
    crit.iterations = 10000;
    solver.setStopCriteria(crit);

    for (int iy1 = 1; iy1 < my; ++iy1) {
    for (int iy2 = iy1; iy2 < my; ++iy2) {
    for (int iz1 = 0; iz1 < mz; ++iz1) {
    for (int iz2 = iz1; iz2 < mz; ++iz2) {

        double scalingy1 = double(iy1) / double(my - 1);
        double y1 = scalingy1 * maxy;
        double scalingy2 = double(iy2) / double(my - 1);
        double y2 = -scalingy2 * maxy;

        double scalingz1 = double(iz1) / double(mz - 1);
        double z1 = scalingz1 * maxz;
        double scalingz2 = double(iz2) / double(mz - 1);
        double z2 = -scalingz2 * maxz;


        Eigen::VectorXd start(6);
        start << 0, y1, z1, 0, y2, z2;
        vec = start;

        electronicWaveFunctionProblem.reset();
        electronicWaveFunctionProblem.putElectronsIntoNuclei(vec, grad);
        solver.minimize(electronicWaveFunctionProblem, vec);
        electronicWaveFunctionProblem.putElectronsIntoNuclei(vec, grad);

        double prec = 0.01;
        if (!vec.isApprox(secondOrderTS1, prec) &&
            !vec.isApprox(secondOrderTS2, prec) &&
            !vec.isApprox(bothRight, prec) &&
            !vec.isApprox(bothLeft, prec) &&
                std::abs(vec[2])<0.7 &&
                std::abs(vec[5])<0.7) {
            std::cout << ParticleCollection(start) << std::endl;
            std::cout << ParticleCollection(vec) << std::endl;

            /*Start Hessian analysis*/
            auto n = ElectronicWaveFunction::getInstance().getNumberOfElectrons() * 3;

            std::cout << "Core ";
            for (auto &it : electronicWaveFunctionProblem.getIndicesOfElectronsAtNuclei()) std::cout << it << " ";
            std::cout << std::endl;
            std::cout << "Free ";
            for (auto &it : electronicWaveFunctionProblem.getIndicesOfElectronsNotAtNuclei())
                std::cout << it << " ";
            std::cout << std::endl;

            Eigen::MatrixXd hess(n, n);
            electronicWaveFunctionProblem.hessian(vec, hess);
            std::cout << hess << std::endl << std::endl;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> selfAdjointEigenSolver(hess, false);

            //auto relevantBlock = hess.block((8 - nsmooth) * 3, (8 - nsmooth) * 3, 3 * nsmooth, 3 * nsmooth);
            //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> selfAdjointEigenSolver(relevantBlock, true);
            //std::cout << relevantBlock << std::endl;

            auto eigenvalues = selfAdjointEigenSolver.eigenvalues();
            std::cout << eigenvalues.transpose() << std::endl;
            //auto eigenvectors = selfAdjointEigenSolver.eigenvectors();
            //std::cout << eigenvectors << std::endl;
            std::cout << std::endl;
            /* End Hessian analysis*/
        }
    }}}}


    if(showGui) {
        QApplication app(argc, argv);
        setlocale(LC_NUMERIC,"C");

        // Prepare the optimization path for visualization
        auto optimizationPath = electronicWaveFunctionProblem.getOptimizationPath();
        ElectronCollections shortenedPath(ElectronCollection(optimizationPath.front(),
                                                             optimizationPath.getSpinTypeCollection()));
        unsigned long nwanted = 300;
        auto skip = 1 + (optimizationPath.length() / nwanted);
        std::cout << "displaying structures with a spacing of " << skip << "." << std::endl;
        for (unsigned long i = 0; i < optimizationPath.length(); i = i + skip) {
            shortenedPath.append(ElectronCollection(optimizationPath.getElectronCollection(i),
                                                    optimizationPath.getSpinTypeCollection()));
        }

        shortenedPath.append(ElectronCollection(vec, optimizationPath.getSpinTypeCollection().asEigenVector()()()));

        // Visualization
        MoleculeWidget moleculeWidget;
        Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

        AtomCollection3D(root, ElectronicWaveFunction::getInstance().getAtomCollection());

        // Plot the starting point
        ElectronCollection3D(root, ElectronCollection(ParticleCollection(vec),
                                                      optimizationPath.getSpinTypeCollection()), false);

        // Plot the optimization path
        ParticleCollectionPath3D(root, shortenedPath);

        return app.exec();
    }
    return 0;
};
