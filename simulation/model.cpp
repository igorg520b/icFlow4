#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;


void icy::Model::Reset() {}

void icy::Model::AssembleAndSolve(SimParams &prms, double timeStep) {}
void icy::Model::GetResultFromSolver(double timeStep) {}
void icy::Model::AcceptTentativeValues(SimParams &prms) {}
void icy::Model::UnsafeUpdateGeometry() {}
void icy::Model::ChangeVisualizationOption(VisOpt option) {}

