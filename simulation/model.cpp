#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"


icy::Model::Model()
{


}


void icy::Model::UnsafeUpdateGeometry()
{
    vtk_update_mutex.lock();
    mesh.UnsafeUpdateGeometry();
    vtk_update_mutex.unlock();
}


void icy::Model::Reset(SimParams &prms)
{
    mesh.Reset(prms.CharacteristicLength);
    UnsafeUpdateGeometry();
}









void icy::Model::InitialGuess(SimParams &prms, double timeStep, double timeStepFactor)
{
    std::size_t nNodes = mesh.allNodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.allNodes[i];
        if(nd->pinned)
        {
            nd->xt = (1/timeStepFactor)*nd->intended_position + (1-1/timeStepFactor)*nd->xn;
            nd->vn=Eigen::Vector2d::Zero();
            nd->x_hat = nd->xn;
        }
        else
        {
            nd->x_hat = nd->xn + timeStep*nd->vn;
            nd->x_hat.y() -= prms.Gravity*timeStep*timeStep;
            nd->xt = nd->xn;
        }
    }
}


bool icy::Model::AssembleAndSolve(SimParams &prms, double timeStep)
{
    qDebug() << "assembling";
    eqOfMotion.ClearAndResize(mesh.freeNodeCount);

    unsigned nElems = mesh.allElems.size();
    unsigned nNodes = mesh.allNodes.size();
    unsigned nInteractions = mesh.collision_interactions.size();

    qDebug() << "assembling - adding elems to sparcity structure " << nElems << " " << mesh.freeNodeCount;

//#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++) {
        qDebug() << "adding elem " << i;
        mesh.allElems[i]->AddToSparsityStructure(eqOfMotion);
    }
/*
    mesh.DetectContactPairs(prms.InteractionDistance);

#pragma omp parallel for
    for(unsigned i=0;i<nInteractions;i++)
        mesh.collision_interactions[i].AddToSparsityStructure(eqOfMotion);
*/
    qDebug() << "assembling - creating structure";
    eqOfMotion.CreateStructure();

    // assemble
    qDebug() << "assembling - elems entries";

    bool mesh_iversion_detected = false;
#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++)
    {
        bool result = mesh.allElems[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);
        if(!result) mesh_iversion_detected = true;
    }

    if(mesh_iversion_detected) return false; // mesh inversion

    qDebug() << "assembling - nodal entries";

#pragma omp parallel for
    for(unsigned i=0;i<nNodes;i++) mesh.allNodes[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);
/*
#pragma omp parallel for
    for(unsigned i=0;i<nInteractions;i++)
        mesh.collision_interactions[i].Evaluate(eqOfMotion, prms, timeStep);

    if(nInteractions==0)
    {
        avgSeparationDistance=-1;
    }
    else
    {
        double distTotal = 0;
#pragma omp parallel for reduction(+:distTotal)
        for(unsigned i=0;i<nInteractions;i++)
            distTotal+=mesh.collision_interactions[i].dist;
        avgSeparationDistance = distTotal/nInteractions;
    }
*/
    // solve
    qDebug() << "solving";
    bool result = eqOfMotion.Solve();

    // pull
    qDebug() << "pulling";
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.allNodes[i];
        Eigen::Vector2d delta_x;
        if(!nd->pinned)
        {
            eqOfMotion.GetTentativeResult(nd->eqId, delta_x);
            nd->xt+=delta_x;
        }
    }
    qDebug() << "pulled";

    return result;
}

void icy::Model::AcceptTentativeValues(double timeStep)
{
    vtk_update_mutex.lock();
    unsigned nNodes = mesh.allNodes.size();
#pragma omp parallel for
    for(unsigned i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.allNodes[i];
        Eigen::Vector2d dx = nd->xt-nd->xn;
        nd->vn = dx/timeStep;
        nd->xn = nd->xt;
    }
    vtk_update_mutex.unlock();
}

