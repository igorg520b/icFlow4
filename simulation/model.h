#ifndef MESHCOLLECTION_H
#define MESHCOLLECTION_H

#include <QFileInfo>
#include <QObject>
#include <QMutex>

#include <vector>
#include <algorithm>
#include <chrono>
#include <unordered_set>




#include "parameters_sim.h"
#include "mesh.h"
#include "equationofmotionsolver.h"

#include <Eigen/Core>

namespace icy { class Model; class Node; class Element;}

class icy::Model : public QObject
{


public:    


    Model();
    void Reset(SimParams &prms);

    void InitialGuess(SimParams &prms, double timeStep, double timeStepFactor);
    bool AssembleAndSolve(SimParams &prms, double timeStep);    // return what solver returns
    void AcceptTentativeValues(double timeStep);
    void UnsafeUpdateGeometry();



    icy::Mesh mesh;



    EquationOfMotionSolver eqOfMotion;
    double avgSeparationDistance = -1;

private:

    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry









signals:
    void requestGeometryUpdate(); // this goes to the main thread, which calls UnsafeUpdateGeometry()
    void propertyChanged();

private:


};

#endif // MESHCOLLECTION_H
