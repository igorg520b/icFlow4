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
#include "equationofmotionsolver.h"

#include <Eigen/Core>

namespace icy { class Model; class Mesh; class Node; class Element; }

class icy::Model : public QObject
{
    Q_OBJECT

public:

    // visualization options
    enum VisOpt { none, elem_area, energy_density, stress_xx, stress_yy, stress_hydrostatic, non_symm_measure,
                ps1, ps2, shear_stress, volume_change};
    Q_ENUM(VisOpt)

    Model();
    ~Model();
    void Reset(SimParams &prms);

    void InitialGuess(SimParams &prms, double timeStep, double timeStepFactor);
    bool AssembleAndSolve(SimParams &prms, double timeStep);    // return what solver returns
    void AcceptTentativeValues(double timeStep);
    void UnsafeUpdateGeometry();
    void ChangeVisualizationOption(icy::Model::VisOpt option);
    void PositionIndenter(double offset);

    icy::Mesh *mesh;

    EquationOfMotionSolver eqOfMotion;
    double avgSeparationDistance = -1;

private:
    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry

signals:
    void requestGeometryUpdate(); // this goes to the main thread, which calls UnsafeUpdateGeometry()
    void propertyChanged();
};

#endif // MESHCOLLECTION_H
