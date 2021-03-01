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

namespace icy { class Model; class Node; class Element;}

class icy::Model : public QObject
{
    Q_OBJECT

public:    

    // visualization options
    enum VisOpt { none, energy_density };
    Q_ENUM(VisOpt)

    void Reset(SimParams &prms);

    void AssembleAndSolve(SimParams &prms, double timeStep);
    void GetResultFromSolver(double timeStep);
    void AcceptTentativeValues(SimParams &prms);
    void UnsafeUpdateGeometry();    // called from the main thread
    void ChangeVisualizationOption(VisOpt option);  // called from the main thread

    icy::Mesh mesh;

private:
    QMutex vtk_update_mutex; // to prevent modifying mesh data while updating VTK representation
    bool vtk_update_requested = false;  // true when signal has been already emitted to update vtk geometry
    VisOpt VisualizingVariable = VisOpt::none;
//    icy::LinearSystem ls;

signals:
    void requestGeometryUpdate(); // this goes to the main thread, which calls UnsafeUpdateGeometry()
    void propertyChanged();
};

#endif // MESHCOLLECTION_H
