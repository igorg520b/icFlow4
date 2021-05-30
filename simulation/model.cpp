#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"
#include "mesh.h"

icy::Model::Model()
{
    mesh = new icy::Mesh();
    Reset(prms);
}

icy::Model::~Model()
{
    delete mesh;
}

void icy::Model::Reset(SimParams &prms)
{
    mesh->Reset(prms.CharacteristicLength);
    UnsafeUpdateGeometry();
}

void icy::Model::Prepare(void)
{
    abortRequested = false;
    timeStepFactor = 1;
}

bool icy::Model::Step(void)
{
    int iter, attempt = 0;
    bool converges=false;
    bool res; // false if matrix is not PSD
    double h;
    do
    {
        h = prms.InitialTimeStep/timeStepFactor; // time step
        InitialGuess(prms, h, timeStepFactor);
        iter = 0;

        do
        {
            if(abortRequested) {Aborting(); return false;}
            res = AssembleAndSolve(prms, h);

            double ratio = iter == 0 ? 0 : eqOfMotion.solution_norm/eqOfMotion.solution_norm_prev;
            converges = (eqOfMotion.solution_norm < prms.ConvergenceCutoff || ratio < prms.ConvergenceEpsilon);

            std::cout << std::scientific << std::setprecision(1);
            std::cout << (currentStep%2 ? ". " : "  ");
            std::cout << attempt << "-";
            std::cout << std::setw(4) <<std::right<< currentStep;
            std::cout << "-"<< std::left << std::setw(2) << iter;
            std::cout << " obj " << std::setw(10) << eqOfMotion.objective_value;
            std::cout << " sln " << std::setw(10) << eqOfMotion.solution_norm;
            if(iter) std::cout << " ra " << std::setw(10) << ratio;
            else std::cout << "tsf " << std::setw(20) << timeStepFactor;
            std::cout << std::endl;
            iter++;
        } while(res && iter < prms.MaxIter && (iter < prms.MinIter || !converges));

        if(!res)
        {
            qDebug() << "Q not PSD";
            attempt++;
            timeStepFactor*=2;
        }
        else if(!converges)
        {
            qDebug() << "sln did not converge";
            attempt++;
            timeStepFactor*=2;
        }
        if(attempt > 20) throw std::runtime_error("could not solve");
    } while (!res || !converges);

    if(timeStepFactor > 1) timeStepFactor /= 1.2;
    if(timeStepFactor < 1) timeStepFactor=1;
    // accept step
    AcceptTentativeValues(h);
    std::cout << std::endl;
    currentStep++;

    emit stepCompleted();
    return(currentStep < prms.MaxSteps);
}

void icy::Model::RequestAbort(void)
{
    abortRequested = true;
}


void icy::Model::Aborting()
{
    //perform any cleanup if step was aborted
    qDebug() << "icy::ModelController::Aborting()";
    abortRequested = false;
    emit stepAborted();
}

void icy::Model::UnsafeUpdateGeometry()
{
    vtk_update_mutex.lock();
    mesh->UnsafeUpdateGeometry();
    vtk_update_mutex.unlock();
}

void icy::Model::ChangeVisualizationOption(icy::Model::VisOpt option)
{
    vtk_update_mutex.lock();
    mesh->ChangeVisualizationOption((int)option);
    vtk_update_mutex.unlock();
}


void icy::Model::PositionIndenter(double offset)
{
//    vtk_update_mutex.lock();
    unsigned n = mesh->indenter.nodes.size();
    Eigen::Vector2d y_direction = Eigen::Vector2d(0,-1.0);
    for(unsigned i=0;i<n;i++)
    {
        icy::Node &nd = mesh->indenter.nodes[i];
        nd.intended_position = nd.x_initial + offset*y_direction;
    }
//    vtk_update_mutex.unlock();
}












void icy::Model::InitialGuess(SimParams &prms, double timeStep, double timeStepFactor)
{
    std::size_t nNodes = mesh->allNodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh->allNodes[i];
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
    eqOfMotion.ClearAndResize(mesh->freeNodeCount);

    unsigned nElems = mesh->allElems.size();
    unsigned nNodes = mesh->allNodes.size();

#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++) mesh->allElems[i]->AddToSparsityStructure(eqOfMotion);

    vtk_update_mutex.lock();
    mesh->DetectContactPairs(prms.InteractionDistance);
    vtk_update_mutex.unlock();
    unsigned nInteractions = mesh->collision_interactions.size();

/*
#pragma omp parallel for
    for(unsigned i=0;i<nInteractions;i++)
        mesh->collision_interactions[i].AddToSparsityStructure(eqOfMotion);
*/
    eqOfMotion.CreateStructure();

    // assemble
    bool mesh_iversion_detected = false;
#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++)
    {
        bool result = mesh->allElems[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);
        if(!result) mesh_iversion_detected = true;
    }

    if(mesh_iversion_detected) return false; // mesh inversion

#pragma omp parallel for
    for(unsigned i=0;i<nNodes;i++) mesh->allNodes[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);
/*
#pragma omp parallel for
    for(unsigned i=0;i<nInteractions;i++)
        mesh->collision_interactions[i].Evaluate(eqOfMotion, prms, timeStep);

    if(nInteractions==0)
    {
        avgSeparationDistance=-1;
    }
    else
    {
        double distTotal = 0;
#pragma omp parallel for reduction(+:distTotal)
        for(unsigned i=0;i<nInteractions;i++)
            distTotal+=mesh->collision_interactions[i].dist;
        avgSeparationDistance = distTotal/nInteractions;
    }
*/
    // solve
    bool result = eqOfMotion.Solve();

    // pull
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh->allNodes[i];
        Eigen::Vector2d delta_x;
        if(!nd->pinned)
        {
            eqOfMotion.GetTentativeResult(nd->eqId, delta_x);
            nd->xt+=delta_x;
        }
    }
    return result;
}

void icy::Model::AcceptTentativeValues(double timeStep)
{
    vtk_update_mutex.lock();
    unsigned nNodes = mesh->allNodes.size();
#pragma omp parallel for
    for(unsigned i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh->allNodes[i];
        Eigen::Vector2d dx = nd->xt-nd->xn;
        nd->vn = dx/timeStep;
        nd->xn = nd->xt;
    }
    vtk_update_mutex.unlock();
}

