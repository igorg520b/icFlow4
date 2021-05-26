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
    timeStepFactor = 1;
}

bool icy::ModelController::Step(void)
{
    int iter, attempt = 0;
    bool converges=false;
    bool res; // false if matrix is not PSD
    double h;
    do
    {
        h = prms.InitialTimeStep/timeStepFactor; // time step
        model.InitialGuess(prms, h, timeStepFactor);
        iter = 0;

        do
        {
            if(abortRequested) {Aborting(); return false;}
            res = model.AssembleAndSolve(prms, h);

            double ratio = iter == 0 ? 0 : model.eqOfMotion.solution_norm/model.eqOfMotion.solution_norm_prev;
            converges = (model.eqOfMotion.solution_norm < prms.ConvergenceCutoff || ratio < prms.ConvergenceEpsilon);

            std::cout << std::scientific << std::setprecision(1);
            std::cout << (currentStep%2 ? ". " : "  ");
            std::cout << attempt << "-";
            std::cout << std::setw(4) <<std::right<< currentStep;
            std::cout << "-"<< std::left << std::setw(2) << iter;
            std::cout << " obj " << std::setw(10) << model.eqOfMotion.objective_value;
            std::cout << " sln " << std::setw(10) << model.eqOfMotion.solution_norm;
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
    model.AcceptTentativeValues(h);
    std::cout << std::endl;
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
