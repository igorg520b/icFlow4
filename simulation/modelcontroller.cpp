#include "modelcontroller.h"

#include <chrono>
#include <thread>
#include <iostream>
#include <iomanip>

icy::ModelController::ModelController()
{
    model.Reset(prms);
}

void icy::ModelController::Prepare(void)
{
    abortRequested = false;
}

bool icy::ModelController::Step(void)
{

    double h = prms.InitialTimeStep; // time step
    model.InitialGuess(prms, h);
    int iter = 0;
    bool converges=false;
    bool res; // false if matrix is not PSD
    do
    {
        if(abortRequested) {Aborting(); return false;}
        res = model.AssembleAndSolve(prms, h);

        double ratio = iter == 0 ? 0 : model.eqOfMotion.solution_norm/model.eqOfMotion.solution_norm_prev;
        converges = (model.eqOfMotion.solution_norm < prms.ConvergenceCutoff || ratio < prms.ConvergenceEpsilon);
        iter++;
        std::cout << std::scientific << std::setprecision(1);
        if(iter>3) std::cout << "* ";
        else if(currentStep%2) std::cout << ". ";
        else std::cout << "  ";
        std::cout << std::setw(4) <<std::right<< currentStep;
        std::cout << "-"<< std::left << std::setw(2) << iter;
        std::cout << " obj " << std::setw(10) << model.eqOfMotion.objective_value;
        std::cout << " sln " << std::setw(10) << model.eqOfMotion.solution_norm;
        if(iter!=1) std::cout << " ra " << std::setw(10) << ratio;
        std::cout << std::endl;

    } while(res && iter < prms.MaxIter && (iter < prms.MinIter || !converges));

    if(!res)
    {
        qDebug() << "Q not PSD";
        throw std::runtime_error("Q not PSD");
    }
    if(!converges)
    {
        qDebug() << "sln did not converge";
        throw std::runtime_error("sln did not converge");
    }

    // accept step
    model.AcceptTentativeValues(h);
    currentStep++;

    emit stepCompleted();
    return(currentStep < prms.MaxSteps);
}

void icy::ModelController::RequestAbort(void)
{
    abortRequested = true;
}


void icy::ModelController::Aborting()
{
    //perform any cleanup if step was aborted
    qDebug() << "icy::ModelController::Aborting()";
    abortRequested = false;
    emit stepAborted();
}
