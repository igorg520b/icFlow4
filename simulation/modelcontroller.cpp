#include "modelcontroller.h"

#include <chrono>
#include <thread>
#include <iostream>

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
    currentStep++;
    std::cout << "ModelControllerTest::Step " << currentStep << std::endl;

    double h = prms.InitialTimeStep; // time step
    model.InitialGuess(prms, h);
    for(int i=0;i<2;i++)
    {
        if(abortRequested) {Aborting(); return false;}
        model.AssembleAndSolve(prms, h);
    }
    model.AcceptTentativeValues(h);
//    std::this_thread::sleep_for(std::chrono::milliseconds(1));    // for testing

    /*
    for(int i=0;i<30;i++) {
        if(!abortRequested) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            std::cout << i << " substep of step " << currentStep << std::endl;
        } else {
            std::cout << "abort requested at " << i << " substep of step " << currentStep << std::endl;
            Aborting();
            break;
        }
    }
    */



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
