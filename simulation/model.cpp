#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"




void icy::Model::Reset(SimParams &prms)
{
    mesh.Reset(prms.CharacteristicLength);
}

void icy::Model::AssembleAndSolve(SimParams &prms, double timeStep) {}
void icy::Model::GetResultFromSolver(double timeStep) {}
void icy::Model::AcceptTentativeValues(SimParams &prms) {}
void icy::Model::UnsafeUpdateGeometry() {}
void icy::Model::ChangeVisualizationOption(VisOpt option) {}

